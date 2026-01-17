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
import json
import argparse
import subprocess
import signal
import threading
from pathlib import Path
from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Optional, List, Tuple
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError

# Default timeout for test operations (in seconds)
TEST_TIMEOUT_SECONDS = 30


class TimeoutError(Exception):
    """Raised when a test times out."""
    pass


def run_with_timeout(func, timeout_seconds=TEST_TIMEOUT_SECONDS):
    """
    Run a function with a timeout. Returns (result, elapsed_ms) or raises TimeoutError.
    Uses ThreadPoolExecutor for cross-platform timeout support.
    """
    start_time = time.perf_counter()
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(func)
        try:
            result = future.result(timeout=timeout_seconds)
            elapsed_ms = (time.perf_counter() - start_time) * 1000
            return result, elapsed_ms
        except FuturesTimeoutError:
            raise TimeoutError(f"Operation timed out after {timeout_seconds} seconds")

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
    genus: int = 0
    stl_file: Optional[str] = None
    error: Optional[str] = None


def compute_genus(m: Manifold) -> int:
    """
    Get the genus of a manifold using Manifold's built-in method.
    Returns -1 for multi-component manifolds.
    """
    if m.is_empty():
        return 0
    return m.genus()


def get_git_info() -> dict:
    """Get current git commit hash and repository URL."""
    info = {"commit_hash": None, "repo_url": None}
    try:
        # Get commit hash
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            info["commit_hash"] = result.stdout.strip()

        # Get remote URL and convert to https
        result = subprocess.run(
            ["git", "remote", "get-url", "origin"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            url = result.stdout.strip()
            # Convert git@github.com:user/repo.git to https://github.com/user/repo
            if url.startswith("git@"):
                url = url.replace(":", "/").replace("git@", "https://")
            if url.endswith(".git"):
                url = url[:-4]
            info["repo_url"] = url
    except Exception:
        pass
    return info


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


def load_stl_as_manifold(filepath: str, merge_tolerance: float = 0.001) -> Optional[Manifold]:
    """
    Load an STL file and return it as a Manifold.
    Merges duplicate vertices within the given tolerance.
    """
    import struct

    try:
        with open(filepath, 'rb') as f:
            header = f.read(80)
            num_triangles = struct.unpack('<I', f.read(4))[0]

            vertices = []
            faces = []

            for i in range(num_triangles):
                struct.unpack('<fff', f.read(12))  # discard normal
                v1 = struct.unpack('<fff', f.read(12))
                v2 = struct.unpack('<fff', f.read(12))
                v3 = struct.unpack('<fff', f.read(12))
                f.read(2)  # attribute byte count

                idx_base = len(vertices)
                vertices.extend([v1, v2, v3])
                faces.append([idx_base, idx_base+1, idx_base+2])

        if not HAVE_NUMPY:
            print("Warning: numpy required for STL loading")
            return None

        verts = np.array(vertices, dtype=np.float32)
        tris = np.array(faces, dtype=np.uint32)

        # Merge duplicate vertices
        unique_verts = []
        vert_map = {}
        for i, v in enumerate(verts):
            key = tuple(np.round(v / merge_tolerance).astype(int))
            if key not in vert_map:
                vert_map[key] = len(unique_verts)
                unique_verts.append(v)

        new_tris = []
        for tri in tris:
            new_tri = []
            for idx in tri:
                key = tuple(np.round(verts[idx] / merge_tolerance).astype(int))
                new_tri.append(vert_map[key])
            new_tris.append(new_tri)

        verts = np.array(unique_verts, dtype=np.float32)
        tris = np.array(new_tris, dtype=np.uint32)

        from manifold3d import Mesh
        mesh = Mesh(vert_properties=verts, tri_verts=tris)
        manifold = Manifold(mesh)

        if manifold.is_empty():
            print(f"Warning: STL loaded but resulted in empty manifold: {filepath}")
            return None

        return manifold
    except Exception as e:
        print(f"Error loading STL {filepath}: {e}")
        return None


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

    # Load spoon STL for Minkowski sum tests (like CGAL's classic example)
    # The spoon is scaled to fit within a reasonable bounding box
    script_dir = Path(__file__).parent
    spoon_path = script_dir / "rice_spoon.stl"
    if spoon_path.exists():
        spoon = load_stl_as_manifold(str(spoon_path))
        if spoon is not None:
            # Scale and center the spoon to fit within a 1.0 unit bounding box
            # Original bounds are roughly -11.5 to 11.5 in x, -40 to 40 in y, 0 to 7 in z
            # Scale to fit in a ~1.0 unit cube, centered
            scale_factor = 1.0 / 80.0  # Scale down by 80x
            spoon = spoon.scale([scale_factor, scale_factor, scale_factor])
            shapes['spoon'] = spoon
    else:
        print(f"Warning: Spoon STL not found at {spoon_path}")

    return shapes


def test_offset_method(shape: Manifold, shape_name: str, delta: float,
                       circular_segments: int, method: OffsetMethod,
                       method_name: str) -> Tuple[TestResult, Optional[Manifold]]:
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
            is_empty=result.is_empty(),
            genus=compute_genus(result)
        ), result
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
        ), None


def test_minkowski_method(shape: Manifold, shape_name: str, delta: float,
                          circular_segments: int) -> Tuple[TestResult, Optional[Manifold]]:
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
            is_empty=result.is_empty(),
            genus=compute_genus(result)
        ), result
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
        ), None


def create_nonconvex_structuring_elements() -> dict:
    """
    Create non-convex structuring elements for NonConvex-NonConvex Minkowski tests.
    These are small non-convex shapes used as the second operand in Minkowski operations.
    """
    elements = {}

    # Small L-shape structuring element
    box1 = Manifold.cube([0.1, 0.03, 0.03], True)
    box2 = Manifold.cube([0.03, 0.1, 0.03], True)
    elements['l_element'] = box1 + box2

    # Small cross/plus shape (non-convex)
    arm_x = Manifold.cube([0.12, 0.03, 0.03], True)
    arm_y = Manifold.cube([0.03, 0.12, 0.03], True)
    arm_z = Manifold.cube([0.03, 0.03, 0.12], True)
    elements['cross_element'] = arm_x + arm_y + arm_z

    # Small notched cube (cube with corner cut out)
    cube = Manifold.cube([0.08, 0.08, 0.08], True)
    corner = Manifold.cube([0.05, 0.05, 0.05], True).translate([0.03, 0.03, 0.03])
    elements['notched_element'] = cube - corner

    # Star element: cube with a point (spike) above each face
    # This matches the classic CGAL Minkowski sum star shape
    half = 0.04  # half-size of the cube
    spike_height = 0.06  # height of each spike from face center
    # Start with the central cube
    star_cube = Manifold.cube([half * 2, half * 2, half * 2], True)
    # Create spikes by hulling each face with a point above it
    # Face vertices for a centered cube: corners at (+/-half, +/-half, +/-half)
    spikes = []
    # +X face spike
    spikes.append(Manifold.hull_points([
        [half, -half, -half], [half, half, -half],
        [half, half, half], [half, -half, half],
        [half + spike_height, 0, 0]
    ]))
    # -X face spike
    spikes.append(Manifold.hull_points([
        [-half, -half, -half], [-half, half, -half],
        [-half, half, half], [-half, -half, half],
        [-half - spike_height, 0, 0]
    ]))
    # +Y face spike
    spikes.append(Manifold.hull_points([
        [-half, half, -half], [half, half, -half],
        [half, half, half], [-half, half, half],
        [0, half + spike_height, 0]
    ]))
    # -Y face spike
    spikes.append(Manifold.hull_points([
        [-half, -half, -half], [half, -half, -half],
        [half, -half, half], [-half, -half, half],
        [0, -half - spike_height, 0]
    ]))
    # +Z face spike
    spikes.append(Manifold.hull_points([
        [-half, -half, half], [half, -half, half],
        [half, half, half], [-half, half, half],
        [0, 0, half + spike_height]
    ]))
    # -Z face spike
    spikes.append(Manifold.hull_points([
        [-half, -half, -half], [half, -half, -half],
        [half, half, -half], [-half, half, -half],
        [0, 0, -half - spike_height]
    ]))
    # Union the cube with all spikes
    star = star_cube
    for spike in spikes:
        star = star + spike
    elements['star_element'] = star

    return elements


def test_nonconvex_minkowski(shape: Manifold, shape_name: str,
                              element: Manifold, element_name: str,
                              is_sum: bool = True,
                              timeout_seconds: int = TEST_TIMEOUT_SECONDS) -> Tuple[TestResult, Optional[Manifold]]:
    """
    Test NonConvex-NonConvex Minkowski sum/difference.
    Both shape and element should be non-convex to exercise the slow path.
    """
    op_name = "sum" if is_sum else "diff"
    test_name = f"{shape_name}_x_{element_name}"
    method_name = f"nc_mink_{op_name}"

    def do_minkowski():
        if is_sum:
            return shape.minkowski_sum(element)
        else:
            return shape.minkowski_difference(element)

    try:
        result, elapsed_ms = run_with_timeout(do_minkowski, timeout_seconds)

        return TestResult(
            name=test_name,
            method=method_name,
            delta=1.0 if is_sum else -1.0,  # Use delta to indicate sum vs diff
            time_ms=elapsed_ms,
            num_verts=result.num_vert(),
            num_tris=result.num_tri(),
            volume=result.volume() if not result.is_empty() else 0.0,
            is_empty=result.is_empty(),
            genus=compute_genus(result)
        ), result
    except TimeoutError as e:
        return TestResult(
            name=test_name,
            method=method_name,
            delta=1.0 if is_sum else -1.0,
            time_ms=timeout_seconds * 1000,
            num_verts=0,
            num_tris=0,
            volume=0.0,
            is_empty=True,
            error=f"TIMEOUT ({timeout_seconds}s)"
        ), None
    except Exception as e:
        return TestResult(
            name=test_name,
            method=method_name,
            delta=1.0 if is_sum else -1.0,
            time_ms=0,
            num_verts=0,
            num_tris=0,
            volume=0.0,
            is_empty=True,
            error=str(e)
        ), None


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
            is_empty=closed.is_empty(),
            genus=compute_genus(closed)
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
            is_empty=opened.is_empty(),
            genus=compute_genus(opened)
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
    print("\n" + "=" * 110)
    print(f"{'Shape':<20} {'Method':<18} {'Delta':>8} {'Time(ms)':>10} "
          f"{'Verts':>8} {'Tris':>8} {'Genus':>6} {'Volume':>12} {'Status':<10}")
    print("=" * 110)

    for r in results:
        status = "ERROR" if r.error else ("EMPTY" if r.is_empty else "OK")
        print(f"{r.name:<20} {r.method:<18} {r.delta:>8.3f} {r.time_ms:>10.1f} "
              f"{r.num_verts:>8} {r.num_tris:>8} {r.genus:>6} {r.volume:>12.4f} {status:<10}")

    print("=" * 110)


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

    print(f"Output directory: {output_path.absolute()}")
    print(f"Circular segments: {circular_segments}")
    print(f"Test deltas: {deltas}")
    print()

    # Create test shapes
    shapes = create_test_shapes()
    print(f"Created {len(shapes)} test shapes: {list(shapes.keys())}")

    # Save original shapes
    print("\nSaving original shapes...", flush=True)
    for name, shape in shapes.items():
        save_stl(shape, str(stl_dir / f"{name}_original.stl"))

    all_results = []
    methods = [
        (OffsetMethod.Minkowski, "minkowski"),
        (OffsetMethod.Simple, "simple"),
        (OffsetMethod.Elegant, "elegant"),
    ]

    # Test each shape with each method and delta
    for shape_name, shape in shapes.items():
        print(f"\nTesting {shape_name}...", flush=True)

        for delta in deltas:
            delta_str = f"d{delta:+.2f}".replace(".", "p").replace("+", "pos").replace("-", "neg")

            for method, method_name in methods:
                print(f"  Starting {method_name} delta={delta:+.2f}...", end=" ", flush=True)
                result, manifold = test_offset_method(shape, shape_name, delta,
                                           circular_segments, method, method_name)

                status = "ERROR" if result.error else ("EMPTY" if result.is_empty else "OK")
                print(f"{result.time_ms:.1f}ms, {result.num_tris} tris, {status}", flush=True)

                # Save results
                stl_filename = f"{shape_name}_{method_name}_{delta_str}.stl"
                if manifold is not None and not result.is_empty and not result.error:
                    try:
                        save_stl(manifold, str(stl_dir / stl_filename))
                        result.stl_file = stl_filename
                    except Exception as e:
                        print(f"    Warning: Could not save {method_name} result: {e}")

                all_results.append(result)

            # Also test direct minkowski for comparison
            print(f"  Starting minkowski_direct delta={delta:+.2f}...", end=" ", flush=True)
            result, _ = test_minkowski_method(shape, shape_name, delta, circular_segments)
            status = "ERROR" if result.error else ("EMPTY" if result.is_empty else "OK")
            print(f"{result.time_ms:.1f}ms, {result.num_tris} tris, {status}", flush=True)
            all_results.append(result)

    # Test morphological operations on fun_shape
    print("\nTesting morphological operations on fun_shape...", flush=True)
    fun_shape = shapes['fun_shape']

    for method, method_name in methods:
        print(f"  Starting morphological {method_name}...", end=" ", flush=True)
        morph_results, final = test_morphological_operations(
            fun_shape, "fun_shape", 0.1, circular_segments, method, method_name)
        print(f"done ({len(morph_results)} results)", flush=True)

        # Set STL filenames for morphological results
        for r in morph_results:
            if "closing" in r.name:
                stl_filename = f"fun_shape_closing_{method_name}.stl"
                if final is not None and not final.is_empty():
                    save_stl(final, str(stl_dir / stl_filename))
                    r.stl_file = stl_filename

        all_results.extend(morph_results)

    # Test NonConvex-NonConvex Minkowski operations
    print("\nTesting NonConvex-NonConvex Minkowski operations...", flush=True)
    nc_elements = create_nonconvex_structuring_elements()

    # Save structuring elements
    for elem_name, elem in nc_elements.items():
        save_stl(elem, str(stl_dir / f"{elem_name}_original.stl"))

    # Select non-convex shapes for testing
    # Note: spoon is excluded because it has 904 triangles, making NC-NC Minkowski
    # operations prohibitively slow (904 * element_tris hull operations)
    nonconvex_shapes = {
        'fun_shape': shapes['fun_shape'],
        'l_shape': shapes['l_shape'],
        # 'hollow_sphere': shapes['hollow_sphere'],  # Also slow due to high tri count
    }

    # Test each non-convex shape with each non-convex structuring element
    for shape_name, shape in nonconvex_shapes.items():
        for elem_name, element in nc_elements.items():
            # Test Minkowski sum (dilation)
            print(f"  Starting {shape_name} + {elem_name} (sum)...", end=" ", flush=True)
            result, manifold = test_nonconvex_minkowski(
                shape, shape_name, element, elem_name, is_sum=True)
            status = "ERROR" if result.error else ("EMPTY" if result.is_empty else "OK")
            print(f"{result.time_ms:.1f}ms, {result.num_tris} tris, {status}", flush=True)

            # Save STL
            stl_filename = f"{shape_name}_x_{elem_name}_sum.stl"
            if manifold is not None and not result.is_empty and not result.error:
                try:
                    save_stl(manifold, str(stl_dir / stl_filename))
                    result.stl_file = stl_filename
                except Exception as e:
                    print(f"    Warning: Could not save result: {e}")
            all_results.append(result)

            # Test Minkowski difference (erosion)
            print(f"  Starting {shape_name} - {elem_name} (diff)...", end=" ", flush=True)
            result, manifold = test_nonconvex_minkowski(
                shape, shape_name, element, elem_name, is_sum=False)
            status = "ERROR" if result.error else ("EMPTY" if result.is_empty else "OK")
            print(f"{result.time_ms:.1f}ms, {result.num_tris} tris, {status}", flush=True)

            # Save STL
            stl_filename = f"{shape_name}_x_{elem_name}_diff.stl"
            if manifold is not None and not result.is_empty and not result.error:
                try:
                    save_stl(manifold, str(stl_dir / stl_filename))
                    result.stl_file = stl_filename
                except Exception as e:
                    print(f"    Warning: Could not save result: {e}")
            all_results.append(result)

    # Print summary table
    print_results_table(all_results)

    # Save results to CSV
    csv_path = output_path / "results.csv"
    with open(csv_path, 'w') as f:
        f.write("shape,method,delta,time_ms,num_verts,num_tris,volume,genus,is_empty,error\n")
        for r in all_results:
            error_str = r.error.replace(',', ';') if r.error else ""
            f.write(f"{r.name},{r.method},{r.delta},{r.time_ms:.2f},"
                   f"{r.num_verts},{r.num_tris},{r.volume:.6f},{r.genus},"
                   f"{r.is_empty},{error_str}\n")
    print(f"\nResults saved to {csv_path}")

    # Generate JSON for web viewer
    json_path = output_path / "comparison_data.json"
    git_info = get_git_info()
    json_data = {
        "timestamp": datetime.now().isoformat(),
        "commit_hash": git_info["commit_hash"],
        "repo_url": git_info["repo_url"],
        "circular_segments": circular_segments,
        "deltas": deltas,
        "results": []
    }
    for r in all_results:
        json_data["results"].append({
            "shape": r.name,
            "method": r.method,
            "delta": r.delta,
            "time_ms": r.time_ms,
            "num_verts": r.num_verts,
            "num_tris": r.num_tris,
            "volume": r.volume,
            "genus": r.genus,
            "is_empty": r.is_empty,
            "stl_file": r.stl_file,
            "error": r.error
        })
    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"JSON data saved to {json_path}")

    # Copy index.html to output directory if not already there
    script_dir = Path(__file__).parent
    index_src = script_dir / "offset_viewer.html"
    index_dst = output_path / "index.html"

    # Generate index.html inline if template doesn't exist
    if not index_dst.exists():
        print(f"Note: Run with web viewer at {output_path}/index.html")

    return all_results


def deploy_to_cloudflare(output_dir: str, project_name: str = "manifold-offset-comparison"):
    """Deploy the output directory to Cloudflare Pages using wrangler."""
    import subprocess
    import shutil

    # Check if wrangler is available
    if not shutil.which("wrangler"):
        print("\nWarning: wrangler not found. Skipping Cloudflare Pages deployment.")
        print("Install with: npm install -g wrangler")
        return False

    print(f"\nDeploying to Cloudflare Pages ({project_name})...")
    try:
        result = subprocess.run(
            ["wrangler", "pages", "deploy", ".",
             f"--project-name={project_name}",
             "--commit-dirty=true",
             "--branch=main"],
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode == 0:
            # Extract URL from output
            for line in result.stdout.split('\n'):
                if 'pages.dev' in line and 'https://' in line:
                    url = line.strip().split()[-1]
                    print(f"Deployed successfully: {url}")
                    break
            else:
                print("Deployed successfully!")
            return True
        else:
            print(f"Deployment failed: {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print("Deployment timed out")
        return False
    except Exception as e:
        print(f"Deployment error: {e}")
        return False


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
    parser.add_argument("--deploy", action="store_true",
                       help="Deploy results to Cloudflare Pages")
    parser.add_argument("--project-name", default="manifold-offset-comparison",
                       help="Cloudflare Pages project name")

    args = parser.parse_args()

    print("=" * 60)
    print("Minkowski vs Offset Comparison Test Suite")
    print("=" * 60)
    print()

    run_tests(args.output, args.deltas, args.segments, not args.quiet)

    if args.deploy:
        deploy_to_cloudflare(args.output, args.project_name)


if __name__ == "__main__":
    main()
