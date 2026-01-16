#!/usr/bin/env python3
"""
Test suite comparing Minkowski sum/difference vs Offset operations for dilation/erosion.

This script evaluates the relative strengths and weaknesses of different approaches:
- Minkowski method: Uses MinkowskiSum/Difference with a sphere (most robust)
- Simple offset: Uses cylinders on convex edges, spheres on convex vertices
- Elegant offset: Uses circular arc wedges on edges, hulled sphere caps on vertices

Based on discussions in PR #669 (https://github.com/elalish/manifold/pull/669)
"""

import os
import sys
import time
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Tuple

# Try to import manifold3d
try:
    import manifold3d
    from manifold3d import Manifold, OffsetMethod
except ImportError:
    print("Error: manifold3d module not found. Please build and install it first.")
    print("Run: pip install -e . from the manifold root directory")
    sys.exit(1)

# Try to import visualization dependencies
try:
    import numpy as np
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False
    print("Warning: numpy not found. Some features may be limited.")

try:
    import trimesh
    HAVE_TRIMESH = True
except Exception as e:
    HAVE_TRIMESH = False
    print(f"Warning: trimesh not available ({type(e).__name__}: {e}). STL export will be limited.")

try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False
    print("Warning: matplotlib not found. Screenshots will not be generated.")


@dataclass
class TestResult:
    """Result of a single offset/minkowski test."""
    name: str
    method: str
    delta: float
    time_ms: float
    num_verts: int
    num_tris: int
    volume: float
    is_empty: bool
    error: Optional[str] = None


def manifold_to_trimesh(m: Manifold) -> Optional['trimesh.Trimesh']:
    """Convert a Manifold to a trimesh object."""
    if not HAVE_TRIMESH or not HAVE_NUMPY:
        return None
    if m.is_empty():
        return None

    mesh = m.to_mesh()
    verts = np.array(mesh.vert_properties[:, :3])
    faces = np.array(mesh.tri_verts).reshape(-1, 3)
    return trimesh.Trimesh(vertices=verts, faces=faces)


def save_stl(m: Manifold, filepath: str) -> bool:
    """Save a Manifold as an STL file."""
    tm = manifold_to_trimesh(m)
    if tm is None:
        return False
    try:
        tm.export(filepath)
        return True
    except Exception as e:
        print(f"Error saving STL: {e}")
        return False


def save_screenshot(m: Manifold, filepath: str, title: str = "", show_wireframe: bool = True) -> bool:
    """Save a screenshot of the manifold with optional wireframe overlay."""
    if not HAVE_MATPLOTLIB or not HAVE_NUMPY:
        return False
    if m.is_empty():
        return False

    try:
        mesh = m.to_mesh()
        verts = np.array(mesh.vert_properties[:, :3])
        faces = np.array(mesh.tri_verts).reshape(-1, 3)

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Create polygon collection for solid faces
        triangles = verts[faces]

        # Add solid faces with slight transparency
        solid_collection = Poly3DCollection(triangles, alpha=0.7)
        solid_collection.set_facecolor('steelblue')
        solid_collection.set_edgecolor('none')
        ax.add_collection3d(solid_collection)

        # Add wireframe overlay
        if show_wireframe:
            wire_collection = Poly3DCollection(triangles, alpha=0.0)
            wire_collection.set_facecolor('none')
            wire_collection.set_edgecolor('darkblue')
            wire_collection.set_linewidth(0.3)
            ax.add_collection3d(wire_collection)

        # Set axis limits
        max_range = np.array([
            verts[:, 0].max() - verts[:, 0].min(),
            verts[:, 1].max() - verts[:, 1].min(),
            verts[:, 2].max() - verts[:, 2].min()
        ]).max() / 2.0

        mid_x = (verts[:, 0].max() + verts[:, 0].min()) * 0.5
        mid_y = (verts[:, 1].max() + verts[:, 1].min()) * 0.5
        mid_z = (verts[:, 2].max() + verts[:, 2].min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        if title:
            ax.set_title(title, fontsize=12)

        # Set viewing angle
        ax.view_init(elev=30, azim=45)

        plt.tight_layout()
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        print(f"Error saving screenshot: {e}")
        return False


def create_test_shapes() -> dict:
    """Create a variety of test shapes for offset comparison."""
    shapes = {}

    # Simple cube
    shapes['cube'] = Manifold.cube([1.0, 1.0, 1.0], True)

    # Sphere
    shapes['sphere'] = Manifold.sphere(0.5, 20)

    # Fun shape from PR #669 discussions: cube with corner cut out
    cube = Manifold.cube([1.0, 1.0, 1.0], True)
    corner_cube = Manifold.cube([1.0, 1.0, 1.0], True).translate([0.5, 0.5, 0.5])
    shapes['fun_shape'] = cube - corner_cube

    # Cylinder
    shapes['cylinder'] = Manifold.cylinder(1.0, 0.4, 0.4, 32, True)

    # Torus-like shape (sphere with smaller sphere subtracted)
    outer = Manifold.sphere(0.5, 20)
    inner = Manifold.sphere(0.3, 16)
    shapes['hollow_sphere'] = outer - inner

    # L-shape
    box1 = Manifold.cube([1.0, 0.3, 0.3], False)
    box2 = Manifold.cube([0.3, 1.0, 0.3], False)
    shapes['l_shape'] = box1 + box2

    # Star-like extrusion (skip if CrossSection not available)
    try:
        from manifold3d import CrossSection
        from math import cos, sin, pi
        points = []
        for i in range(10):
            angle = i * 2 * pi / 10
            r = 0.5 if i % 2 == 0 else 0.25
            points.append([r * cos(angle), r * sin(angle)])
        star_cs = CrossSection([points])
        shapes['star'] = Manifold.extrude(star_cs, 0.3)
    except Exception as e:
        print(f"Warning: Could not create star shape: {e}")

    return shapes


def test_offset_method(shape: Manifold, shape_name: str, delta: float,
                       circular_segments: int, method: OffsetMethod,
                       method_name: str) -> TestResult:
    """Test a single offset operation and return results."""
    try:
        start_time = time.perf_counter()
        result = shape.offset(delta, circular_segments, method)
        elapsed_ms = (time.perf_counter() - start_time) * 1000

        return TestResult(
            name=shape_name,
            method=method_name,
            delta=delta,
            time_ms=elapsed_ms,
            num_verts=result.num_vert(),
            num_tris=result.num_tri(),
            volume=result.volume() if not result.is_empty() else 0.0,
            is_empty=result.is_empty()
        )
    except Exception as e:
        return TestResult(
            name=shape_name,
            method=method_name,
            delta=delta,
            time_ms=0,
            num_verts=0,
            num_tris=0,
            volume=0.0,
            is_empty=True,
            error=str(e)
        )


def test_minkowski_method(shape: Manifold, shape_name: str, delta: float,
                          circular_segments: int) -> TestResult:
    """Test Minkowski sum/difference directly for comparison."""
    try:
        sphere = Manifold.sphere(abs(delta), circular_segments)

        start_time = time.perf_counter()
        if delta > 0:
            result = shape.minkowski_sum(sphere)
        else:
            result = shape.minkowski_difference(sphere)
        elapsed_ms = (time.perf_counter() - start_time) * 1000

        return TestResult(
            name=shape_name,
            method="minkowski_direct",
            delta=delta,
            time_ms=elapsed_ms,
            num_verts=result.num_vert(),
            num_tris=result.num_tri(),
            volume=result.volume() if not result.is_empty() else 0.0,
            is_empty=result.is_empty()
        )
    except Exception as e:
        return TestResult(
            name=shape_name,
            method="minkowski_direct",
            delta=delta,
            time_ms=0,
            num_verts=0,
            num_tris=0,
            volume=0.0,
            is_empty=True,
            error=str(e)
        )


def test_morphological_operations(shape: Manifold, shape_name: str, delta: float,
                                  circular_segments: int, method: OffsetMethod,
                                  method_name: str) -> Tuple[List[TestResult], Optional[Manifold]]:
    """
    Test morphological opening/closing operations.

    Opening: erosion followed by dilation (removes small protrusions)
    Closing: dilation followed by erosion (fills small holes)
    """
    results = []

    # Morphological closing: dilate then erode
    try:
        start_time = time.perf_counter()
        dilated = shape.offset(delta, circular_segments, method)
        closed = dilated.offset(-delta, circular_segments, method)
        elapsed_ms = (time.perf_counter() - start_time) * 1000

        results.append(TestResult(
            name=f"{shape_name}_closing",
            method=method_name,
            delta=delta,
            time_ms=elapsed_ms,
            num_verts=closed.num_vert(),
            num_tris=closed.num_tri(),
            volume=closed.volume() if not closed.is_empty() else 0.0,
            is_empty=closed.is_empty()
        ))
        final_result = closed
    except Exception as e:
        results.append(TestResult(
            name=f"{shape_name}_closing",
            method=method_name,
            delta=delta,
            time_ms=0, num_verts=0, num_tris=0, volume=0.0, is_empty=True,
            error=str(e)
        ))
        final_result = None

    # Morphological opening: erode then dilate
    try:
        start_time = time.perf_counter()
        eroded = shape.offset(-delta, circular_segments, method)
        opened = eroded.offset(delta, circular_segments, method)
        elapsed_ms = (time.perf_counter() - start_time) * 1000

        results.append(TestResult(
            name=f"{shape_name}_opening",
            method=method_name,
            delta=delta,
            time_ms=elapsed_ms,
            num_verts=opened.num_vert(),
            num_tris=opened.num_tri(),
            volume=opened.volume() if not opened.is_empty() else 0.0,
            is_empty=opened.is_empty()
        ))
    except Exception as e:
        results.append(TestResult(
            name=f"{shape_name}_opening",
            method=method_name,
            delta=delta,
            time_ms=0, num_verts=0, num_tris=0, volume=0.0, is_empty=True,
            error=str(e)
        ))

    return results, final_result


def print_results_table(results: List[TestResult]):
    """Print results in a formatted table."""
    print("\n" + "=" * 100)
    print(f"{'Shape':<20} {'Method':<18} {'Delta':>8} {'Time(ms)':>10} "
          f"{'Verts':>8} {'Tris':>8} {'Volume':>12} {'Status':<10}")
    print("=" * 100)

    for r in results:
        status = "ERROR" if r.error else ("EMPTY" if r.is_empty else "OK")
        print(f"{r.name:<20} {r.method:<18} {r.delta:>8.3f} {r.time_ms:>10.1f} "
              f"{r.num_verts:>8} {r.num_tris:>8} {r.volume:>12.4f} {status:<10}")

    print("=" * 100)


def run_tests(output_dir: str, deltas: List[float] = None,
              circular_segments: int = 25, verbose: bool = True):
    """Run the full test suite."""
    if deltas is None:
        deltas = [0.05, 0.1, -0.05, -0.1]

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    stl_dir = output_path / "stl"
    stl_dir.mkdir(exist_ok=True)
    img_dir = output_path / "screenshots"
    img_dir.mkdir(exist_ok=True)

    print(f"Output directory: {output_path.absolute()}")
    print(f"Circular segments: {circular_segments}")
    print(f"Test deltas: {deltas}")
    print()

    # Create test shapes
    shapes = create_test_shapes()
    print(f"Created {len(shapes)} test shapes: {list(shapes.keys())}")

    # Save original shapes
    print("\nSaving original shapes...")
    for name, shape in shapes.items():
        save_stl(shape, str(stl_dir / f"{name}_original.stl"))
        save_screenshot(shape, str(img_dir / f"{name}_original.png"),
                       f"{name} (original)", show_wireframe=True)

    all_results = []
    methods = [
        (OffsetMethod.Minkowski, "minkowski"),
        (OffsetMethod.Simple, "simple"),
        (OffsetMethod.Elegant, "elegant"),
    ]

    # Test each shape with each method and delta
    for shape_name, shape in shapes.items():
        print(f"\nTesting {shape_name}...")

        for delta in deltas:
            delta_str = f"d{delta:+.2f}".replace(".", "p").replace("+", "pos").replace("-", "neg")

            for method, method_name in methods:
                result = test_offset_method(shape, shape_name, delta,
                                           circular_segments, method, method_name)
                all_results.append(result)

                if verbose:
                    status = "ERROR" if result.error else ("EMPTY" if result.is_empty else "OK")
                    print(f"  {method_name} delta={delta:+.2f}: {result.time_ms:.1f}ms, "
                          f"{result.num_tris} tris, {status}")

                # Save results
                if not result.is_empty and not result.error:
                    try:
                        offset_result = shape.offset(delta, circular_segments, method)
                        save_stl(offset_result,
                                str(stl_dir / f"{shape_name}_{method_name}_{delta_str}.stl"))
                        save_screenshot(offset_result,
                                       str(img_dir / f"{shape_name}_{method_name}_{delta_str}.png"),
                                       f"{shape_name} {method_name} delta={delta:+.2f}",
                                       show_wireframe=True)
                    except Exception as e:
                        print(f"    Warning: Could not save {method_name} result: {e}")

            # Also test direct minkowski for comparison
            result = test_minkowski_method(shape, shape_name, delta, circular_segments)
            all_results.append(result)

    # Test morphological operations on fun_shape
    print("\nTesting morphological operations on fun_shape...")
    fun_shape = shapes['fun_shape']

    for method, method_name in methods:
        morph_results, final = test_morphological_operations(
            fun_shape, "fun_shape", 0.1, circular_segments, method, method_name)
        all_results.extend(morph_results)

        if final is not None and not final.is_empty():
            save_stl(final, str(stl_dir / f"fun_shape_closing_{method_name}.stl"))
            save_screenshot(final,
                           str(img_dir / f"fun_shape_closing_{method_name}.png"),
                           f"fun_shape closing {method_name}",
                           show_wireframe=True)

    # Print summary table
    print_results_table(all_results)

    # Save results to CSV
    csv_path = output_path / "results.csv"
    with open(csv_path, 'w') as f:
        f.write("shape,method,delta,time_ms,num_verts,num_tris,volume,is_empty,error\n")
        for r in all_results:
            error_str = r.error.replace(',', ';') if r.error else ""
            f.write(f"{r.name},{r.method},{r.delta},{r.time_ms:.2f},"
                   f"{r.num_verts},{r.num_tris},{r.volume:.6f},"
                   f"{r.is_empty},{error_str}\n")
    print(f"\nResults saved to {csv_path}")

    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Compare Minkowski and Offset approaches for dilation/erosion")
    parser.add_argument("-o", "--output", default="offset_test_results",
                       help="Output directory for STLs and screenshots")
    parser.add_argument("-s", "--segments", type=int, default=25,
                       help="Circular segments for sphere approximation")
    parser.add_argument("-d", "--deltas", type=float, nargs="+",
                       default=[0.05, 0.1, -0.05, -0.1],
                       help="Offset delta values to test")
    parser.add_argument("-q", "--quiet", action="store_true",
                       help="Reduce output verbosity")

    args = parser.parse_args()

    print("=" * 60)
    print("Minkowski vs Offset Comparison Test Suite")
    print("=" * 60)
    print()

    run_tests(args.output, args.deltas, args.segments, not args.quiet)


if __name__ == "__main__":
    main()
