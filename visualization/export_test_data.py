#!/usr/bin/env python3
"""
Export tetrahedralization test results to JSON for Three.js visualization.
"""

import json
import numpy as np
import sys
import pathlib
import importlib.util
from time import time


def load_manifold3d():
    """Load manifold3d, preferring local build over installed version."""
    script_dir = pathlib.Path(__file__).parent
    potential_paths = [
        script_dir / ".." / "build-tet" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / "build" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / "build-debug" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
        script_dir / ".." / "build-release" / "bindings" / "python" / "manifold3d.cpython-39-x86_64-linux-gnu.so",
    ]

    for path in potential_paths:
        if path.exists():
            spec = importlib.util.spec_from_file_location('manifold3d', str(path.resolve()))
            if spec:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                return module

    import manifold3d
    return manifold3d


def compute_tet_quality(v0, v1, v2, v3):
    """Compute tetrahedron quality metric (0-1)."""
    d0 = v1 - v0
    d1 = v2 - v0
    d2 = v3 - v0
    d3 = v2 - v1
    d4 = v3 - v2
    d5 = v1 - v3

    s0, s1, s2 = np.linalg.norm(d0), np.linalg.norm(d1), np.linalg.norm(d2)
    s3, s4, s5 = np.linalg.norm(d3), np.linalg.norm(d4), np.linalg.norm(d5)

    ms = (s0**2 + s1**2 + s2**2 + s3**2 + s4**2 + s5**2) / 6.0
    rms = np.sqrt(ms)

    if rms == 0:
        return 0.0

    s = 12.0 / np.sqrt(2.0)
    vol = np.dot(d0, np.cross(d1, d2)) / 6.0
    return float(np.clip(s * vol / (rms ** 3), 0, 1))


def tet_mesh_to_dict(tet_mesh, name, description, surface_mesh=None):
    """Convert TetMesh to JSON-serializable dict."""
    vertices = tet_mesh.vert_pos.tolist()
    tetrahedra = tet_mesh.tet_verts.tolist()

    # Compute qualities
    qualities = []
    for tet in tetrahedra:
        v0 = np.array(vertices[tet[0]])
        v1 = np.array(vertices[tet[1]])
        v2 = np.array(vertices[tet[2]])
        v3 = np.array(vertices[tet[3]])
        qualities.append(compute_tet_quality(v0, v1, v2, v3))

    result = {
        "name": name,
        "description": description,
        "vertices": vertices,
        "tetrahedra": tetrahedra,
        "qualities": qualities,
        "stats": {
            "numVertices": len(vertices),
            "numTetrahedra": len(tetrahedra),
            "minQuality": float(min(qualities)) if qualities else 0,
            "maxQuality": float(max(qualities)) if qualities else 0,
            "meanQuality": float(np.mean(qualities)) if qualities else 0,
        }
    }

    # Add surface mesh if provided
    if surface_mesh is not None:
        mesh = surface_mesh.to_mesh()
        result["surface"] = {
            "vertices": mesh.vert_properties[:, :3].tolist(),
            "triangles": mesh.tri_verts.tolist()
        }

    return result


def generate_test_cases(m3d):
    """Generate all test cases and return as list of dicts."""
    test_cases = []

    # 1. Random points (unconstrained)
    print("Generating: Random Points...")
    np.random.seed(42)
    points = [(float(x), float(y), float(z)) for x, y, z in np.random.randn(50, 3)]
    t0 = time()
    tet_mesh = m3d.delaunay_tetrahedralization(points)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Random Points", f"50 random 3D points, unconstrained Delaunay ({elapsed*1000:.1f}ms)"),
        "type": "unconstrained"
    })

    # 2. Cube
    print("Generating: Cube...")
    cube = m3d.Manifold.cube((2, 2, 2), center=True)
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(cube)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Cube", f"Unit cube, constrained ({elapsed*1000:.1f}ms)", cube),
        "type": "constrained"
    })

    # 3. Sphere
    print("Generating: Sphere...")
    sphere = m3d.Manifold.sphere(1.0, 16)
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(sphere)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Sphere", f"Geodesic sphere (16 segments), constrained ({elapsed*1000:.1f}ms)", sphere),
        "type": "constrained"
    })

    # 4. Cylinder
    print("Generating: Cylinder...")
    cylinder = m3d.Manifold.cylinder(2.0, 0.5, circular_segments=16)
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(cylinder)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Cylinder", f"Cylinder (height=2, radius=0.5), constrained ({elapsed*1000:.1f}ms)", cylinder),
        "type": "constrained"
    })

    # 5. Torus
    print("Generating: Torus...")
    m3d.set_circular_segments(16)
    circle = m3d.CrossSection.circle(0.3)
    circle = circle.translate((1.0, 0.0))
    torus = circle.revolve(circular_segments=16)
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(torus)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Torus", f"Torus (major=1.0, minor=0.3), constrained ({elapsed*1000:.1f}ms)", torus),
        "type": "constrained"
    })

    # 6. Boolean difference (cube - sphere)
    print("Generating: Boolean Difference...")
    cube = m3d.Manifold.cube((2, 2, 2), center=True)
    sphere = m3d.Manifold.sphere(1.3, 16)
    bool_result = cube - sphere
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(bool_result)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Boolean Difference", f"Cube minus sphere, constrained ({elapsed*1000:.1f}ms)", bool_result),
        "type": "constrained"
    })

    # 7. Extruded L-shape
    print("Generating: Extruded L-Shape...")
    square1 = m3d.CrossSection.square((2, 1))
    square2 = m3d.CrossSection.square((1, 2))
    l_shape = square1 + square2
    extruded = l_shape.extrude(1.5)
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(extruded)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Extruded L-Shape", f"L-shaped extrusion, constrained ({elapsed*1000:.1f}ms)", extruded),
        "type": "constrained"
    })

    # 8. High-genus shape (cube with holes)
    print("Generating: High-Genus Shape...")
    cube = m3d.Manifold.cube((3, 3, 3), center=True)
    hole_x = m3d.Manifold.cylinder(4, 0.3, circular_segments=8).rotate(y_degrees=90)
    hole_y = m3d.Manifold.cylinder(4, 0.3, circular_segments=8).rotate(x_degrees=90)
    hole_z = m3d.Manifold.cylinder(4, 0.3, circular_segments=8).translate((0, 0, -2))
    high_genus = cube - hole_x - hole_y - hole_z
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(high_genus)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "High-Genus Shape", f"Cube with 3 cylindrical holes (genus={high_genus.genus()}), constrained ({elapsed*1000:.1f}ms)", high_genus),
        "type": "constrained"
    })

    # 9. Tetrahedron
    print("Generating: Tetrahedron...")
    tetrahedron = m3d.Manifold.tetrahedron()
    t0 = time()
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(tetrahedron)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Tetrahedron", f"Single tetrahedron, constrained ({elapsed*1000:.1f}ms)", tetrahedron),
        "type": "constrained"
    })

    # 10. Large random point cloud
    print("Generating: Large Point Cloud...")
    np.random.seed(123)
    points = [(float(x), float(y), float(z)) for x, y, z in np.random.randn(200, 3)]
    t0 = time()
    tet_mesh = m3d.delaunay_tetrahedralization(points)
    elapsed = time() - t0
    test_cases.append({
        **tet_mesh_to_dict(tet_mesh, "Large Point Cloud", f"200 random 3D points, unconstrained ({elapsed*1000:.1f}ms)"),
        "type": "unconstrained"
    })

    return test_cases


def main():
    print("Loading manifold3d...")
    m3d = load_manifold3d()

    if not hasattr(m3d, 'delaunay_tetrahedralization'):
        print("Error: manifold3d does not have tetrahedralization support")
        sys.exit(1)

    print(f"Loaded from: {m3d.__file__}")

    print("\nGenerating test cases...")
    test_cases = generate_test_cases(m3d)

    output_path = pathlib.Path(__file__).parent / "public" / "test_data.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nWriting to {output_path}...")
    with open(output_path, 'w') as f:
        json.dump({
            "generated": time(),
            "testCases": test_cases
        }, f)

    # Calculate file size
    file_size = output_path.stat().st_size
    print(f"Exported {len(test_cases)} test cases ({file_size / 1024:.1f} KB)")

    # Print summary
    print("\nTest Cases Summary:")
    print("-" * 60)
    for tc in test_cases:
        stats = tc['stats']
        print(f"  {tc['name']:25s} | {stats['numTetrahedra']:5d} tets | "
              f"quality: {stats['minQuality']:.2f}-{stats['maxQuality']:.2f} (mean: {stats['meanQuality']:.2f})")


if __name__ == "__main__":
    main()
