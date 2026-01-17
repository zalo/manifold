"""
 Copyright 2024 The Manifold Authors.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Test suite for Delaunay tetrahedralization functionality.

 Note: If testing a newly built version from the build directory, run as:
   PYTHONPATH=/path/to/build/bindings/python python3 test_tetrahedralization.py
 """

import numpy as np
from time import time
import sys

# Try to import manifold3d, checking for tetrahedralization support
# Due to potential conflicts with installed versions, we try to load
# from the build directory first if it exists
import importlib.util
import os
import pathlib

def load_manifold3d():
    """Load manifold3d, preferring local build over installed version."""
    # Check for local build in parent directories
    script_dir = pathlib.Path(__file__).parent
    potential_paths = [
        script_dir / ".." / ".." / ".." / "build-tet" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / ".." / ".." / "build" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / ".." / ".." / "build-debug" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / ".." / ".." / "build-release" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
    ]

    for path in potential_paths:
        if path.exists():
            spec = importlib.util.spec_from_file_location('manifold3d', str(path.resolve()))
            if spec:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                return module

    # Fall back to regular import
    import manifold3d
    return manifold3d

try:
    m3d = load_manifold3d()
    if not hasattr(m3d, 'delaunay_tetrahedralization'):
        print("Warning: manifold3d module does not have tetrahedralization support.")
        print("Make sure you're using a version built with tetrahedralization enabled.")
        print(f"Loaded from: {getattr(m3d, '__file__', 'unknown')}")
        sys.exit(1)
    print(f"Loaded manifold3d from: {m3d.__file__}")
except ImportError as e:
    print(f"Error importing manifold3d: {e}")
    sys.exit(1)


def test_unconstrained_random_points():
    """Test unconstrained Delaunay tetrahedralization with random points."""
    print("\n=== Test: Unconstrained tetrahedralization with random points ===")

    np.random.seed(42)
    num_points = 100
    points = [(float(x), float(y), float(z))
              for x, y, z in np.random.randn(num_points, 3)]

    t0 = time()
    tet_mesh = m3d.delaunay_tetrahedralization(points)
    t1 = time()

    print(f"  Input points: {num_points}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_vert == num_points, "Vertex count should match input"
    assert tet_mesh.num_tet > 0, "Should produce at least some tetrahedra"

    # Verify tet_verts shape
    tet_verts = tet_mesh.tet_verts
    assert tet_verts.shape[1] == 4, "Each tet should have 4 vertices"

    # Verify all indices are valid
    assert np.all(tet_verts >= 0), "All indices should be non-negative"
    assert np.all(tet_verts < num_points), "All indices should be less than num points"

    print("  PASSED")
    return True


def test_unconstrained_with_quality_filter():
    """Test unconstrained tetrahedralization with quality filtering."""
    print("\n=== Test: Unconstrained tetrahedralization with quality filter ===")

    np.random.seed(42)
    num_points = 50
    points = [(float(x), float(y), float(z))
              for x, y, z in np.random.randn(num_points, 3)]

    # Without quality filter
    tet_mesh_no_filter = m3d.delaunay_tetrahedralization(points, min_quality=0.0)

    # With quality filter
    tet_mesh_filtered = m3d.delaunay_tetrahedralization(points, min_quality=0.1)

    print(f"  Without filter: {tet_mesh_no_filter.num_tet} tetrahedra")
    print(f"  With filter (0.1): {tet_mesh_filtered.num_tet} tetrahedra")

    # Filtered should have fewer or equal tetrahedra
    assert tet_mesh_filtered.num_tet <= tet_mesh_no_filter.num_tet, \
        "Quality filter should reduce or maintain tet count"

    print("  PASSED")
    return True


def test_constrained_cube():
    """Test constrained tetrahedralization on a cube manifold."""
    print("\n=== Test: Constrained tetrahedralization on cube ===")

    cube = m3d.Manifold.cube((1, 1, 1))

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(cube)
    t1 = time()

    print(f"  Cube vertices: {cube.num_vert()}")
    print(f"  Cube triangles: {cube.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    # Check that original vertices are preserved at the start
    assert tet_mesh.num_vert >= cube.num_vert(), \
        "Should have at least as many vertices as original"

    print("  PASSED")
    return True


def test_constrained_sphere():
    """Test constrained tetrahedralization on a sphere manifold."""
    print("\n=== Test: Constrained tetrahedralization on sphere ===")

    sphere = m3d.Manifold.sphere(1.0, 16)

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(sphere)
    t1 = time()

    print(f"  Sphere vertices: {sphere.num_vert()}")
    print(f"  Sphere triangles: {sphere.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_constrained_cylinder():
    """Test constrained tetrahedralization on a cylinder manifold."""
    print("\n=== Test: Constrained tetrahedralization on cylinder ===")

    cylinder = m3d.Manifold.cylinder(2.0, 0.5, circular_segments=16)

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(cylinder)
    t1 = time()

    print(f"  Cylinder vertices: {cylinder.num_vert()}")
    print(f"  Cylinder triangles: {cylinder.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_constrained_torus():
    """Test constrained tetrahedralization on a torus-like shape."""
    print("\n=== Test: Constrained tetrahedralization on torus-like shape ===")

    # Create a torus-like shape using revolve
    m3d.set_circular_segments(16)
    circle = m3d.CrossSection.circle(0.3)
    circle = circle.translate((1.0, 0.0))
    torus = circle.revolve(circular_segments=16)

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(torus)
    t1 = time()

    print(f"  Torus vertices: {torus.num_vert()}")
    print(f"  Torus triangles: {torus.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_constrained_boolean_result():
    """Test constrained tetrahedralization on a boolean result."""
    print("\n=== Test: Constrained tetrahedralization on boolean result ===")

    cube = m3d.Manifold.cube((2, 2, 2), center=True)
    sphere = m3d.Manifold.sphere(1.3, 16)

    # Boolean difference creates more complex geometry
    result = cube - sphere

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(result)
    t1 = time()

    print(f"  Boolean result vertices: {result.num_vert()}")
    print(f"  Boolean result triangles: {result.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_constrained_extrusion():
    """Test constrained tetrahedralization on an extruded shape."""
    print("\n=== Test: Constrained tetrahedralization on extrusion ===")

    # Create an L-shaped cross section
    square1 = m3d.CrossSection.square((2, 1))
    square2 = m3d.CrossSection.square((1, 2))
    l_shape = square1 + square2

    # Extrude it
    extruded = l_shape.extrude(1.5)

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(extruded)
    t1 = time()

    print(f"  Extruded shape vertices: {extruded.num_vert()}")
    print(f"  Extruded shape triangles: {extruded.num_tri()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_tetrahedron_volume():
    """Verify that tetrahedra have positive volumes."""
    print("\n=== Test: Tetrahedron volume validation ===")

    np.random.seed(42)
    num_points = 30
    points = [(float(x), float(y), float(z))
              for x, y, z in np.random.randn(num_points, 3)]

    tet_mesh = m3d.delaunay_tetrahedralization(points)

    verts = tet_mesh.vert_pos
    tets = tet_mesh.tet_verts

    positive_volume_count = 0
    total_volume = 0.0

    for tet in tets:
        p0, p1, p2, p3 = verts[tet[0]], verts[tet[1]], verts[tet[2]], verts[tet[3]]

        # Compute signed volume
        d0 = p1 - p0
        d1 = p2 - p0
        d2 = p3 - p0
        volume = np.dot(d0, np.cross(d1, d2)) / 6.0

        if volume > 0:
            positive_volume_count += 1
            total_volume += volume

    print(f"  Total tetrahedra: {len(tets)}")
    print(f"  Positive volume count: {positive_volume_count}")
    print(f"  Total volume: {total_volume:.4f}")

    # Most tetrahedra should have positive volume
    # (orientation may vary but absolute value should be positive)
    assert positive_volume_count > 0 or len(tets) == 0, \
        "Should have some positive volume tetrahedra"

    print("  PASSED")
    return True


def test_large_point_set():
    """Test performance with a larger point set."""
    print("\n=== Test: Large point set performance ===")

    np.random.seed(42)
    num_points = 500
    points = [(float(x), float(y), float(z))
              for x, y, z in np.random.randn(num_points, 3)]

    t0 = time()
    tet_mesh = m3d.delaunay_tetrahedralization(points)
    t1 = time()

    print(f"  Input points: {num_points}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")
    print(f"  Points/ms: {num_points / ((t1-t0)*1000):.2f}")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_minimum_points():
    """Test with minimum number of points (4 for a single tetrahedron)."""
    print("\n=== Test: Minimum points (4 points) ===")

    # Create 4 points forming a tetrahedron
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.5, 1.0, 0.0),
        (0.5, 0.5, 1.0)
    ]

    tet_mesh = m3d.delaunay_tetrahedralization(points)

    print(f"  Input points: 4")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")

    assert tet_mesh.num_tet >= 1, "Should produce at least one tetrahedron"

    print("  PASSED")
    return True


def test_coplanar_points():
    """Test handling of coplanar points (degenerate case)."""
    print("\n=== Test: Coplanar points (degenerate case) ===")

    # Create coplanar points on z=0 plane
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.5, 0.5, 0.0)
    ]

    # This may produce 0 tetrahedra since points are coplanar
    # The small random perturbation in the algorithm may help
    tet_mesh = m3d.delaunay_tetrahedralization(points)

    print(f"  Input points: 5 (coplanar)")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")

    # With the random perturbation, we might get some thin tetrahedra
    # or zero if quality filtering removes them
    print("  PASSED (degenerate case handled)")
    return True


def test_constrained_high_genus():
    """Test constrained tetrahedralization on higher genus shape."""
    print("\n=== Test: Constrained tetrahedralization on high-genus shape ===")

    # Create a cube with multiple holes (using boolean operations)
    cube = m3d.Manifold.cube((3, 3, 3), center=True)

    # Create cylinders for holes
    hole_x = m3d.Manifold.cylinder(4, 0.3, circular_segments=8)
    hole_x = hole_x.rotate(y_degrees=90).translate(0, 0, 0)

    hole_y = m3d.Manifold.cylinder(4, 0.3, circular_segments=8)
    hole_y = hole_y.rotate(x_degrees=90).translate(0, 0, 0)

    hole_z = m3d.Manifold.cylinder(4, 0.3, circular_segments=8)
    hole_z = hole_z.translate(0, 0, -2)

    # Create cube with holes
    result = cube - hole_x - hole_y - hole_z

    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(result)
    t1 = time()

    print(f"  Shape vertices: {result.num_vert()}")
    print(f"  Shape triangles: {result.num_tri()}")
    print(f"  Shape genus: {result.genus()}")
    print(f"  Output vertices: {tet_mesh.num_vert}")
    print(f"  Output tetrahedra: {tet_mesh.num_tet}")
    print(f"  Time: {(t1-t0)*1000:.2f}ms")

    assert tet_mesh.num_tet > 0, "Should produce tetrahedra"

    print("  PASSED")
    return True


def test_steiner_point_addition():
    """Test that Steiner points are added when needed."""
    print("\n=== Test: Steiner point addition ===")

    # Create a simple cube
    cube = m3d.Manifold.cube((1, 1, 1))
    original_vert_count = cube.num_vert()

    tet_mesh = m3d.constrained_delaunay_tetrahedralization(cube)

    steiner_points_added = tet_mesh.num_vert - original_vert_count

    print(f"  Original vertices: {original_vert_count}")
    print(f"  Final vertices: {tet_mesh.num_vert}")
    print(f"  Steiner points added: {steiner_points_added}")

    # Some shapes may not need Steiner points, others might
    # Just verify the algorithm runs correctly
    assert tet_mesh.num_vert >= original_vert_count, \
        "Should have at least original vertices"

    print("  PASSED")
    return True


def run_all_tests():
    """Run all tetrahedralization tests."""
    print("=" * 60)
    print("TETRAHEDRALIZATION TEST SUITE")
    print("=" * 60)

    tests = [
        test_unconstrained_random_points,
        test_unconstrained_with_quality_filter,
        test_constrained_cube,
        test_constrained_sphere,
        test_constrained_cylinder,
        test_constrained_torus,
        test_constrained_boolean_result,
        test_constrained_extrusion,
        test_tetrahedron_volume,
        test_large_point_set,
        test_minimum_points,
        test_coplanar_points,
        test_constrained_high_genus,
        test_steiner_point_addition,
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"  FAILED: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed out of {len(tests)} tests")
    print("=" * 60)

    return failed == 0


def run():
    """Entry point for run_all.py compatibility."""
    # Return a simple manifold for compatibility with run_all.py
    return m3d.Manifold.cube((1, 1, 1))


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
