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

 Visualization utilities for Delaunay tetrahedralization results.

 Requires: trimesh, numpy
 Optional: pyglet (for interactive viewing)

 Usage:
   python visualize_tetrahedralization.py
 """

import numpy as np
import sys
import importlib.util
import pathlib


def load_manifold3d():
    """Load manifold3d, preferring local build over installed version."""
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

    import manifold3d
    return manifold3d


try:
    import trimesh
except ImportError:
    print("Error: trimesh is required for visualization.")
    print("Install with: pip install trimesh")
    sys.exit(1)

m3d = load_manifold3d()
if not hasattr(m3d, 'delaunay_tetrahedralization'):
    print("Error: manifold3d does not have tetrahedralization support.")
    sys.exit(1)


def compute_tet_quality(v0, v1, v2, v3):
    """
    Compute tetrahedron quality metric.
    Returns a value between 0 (degenerate) and 1 (regular tetrahedron).
    """
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
    return s * vol / (rms ** 3)


def visualize_tetrahedra(tet_mesh, scale=0.8, color_mode='random',
                         min_volume=1e-6, show_edges=False):
    """
    Visualize tetrahedra using trimesh with exploded view.

    Parameters:
    -----------
    tet_mesh : TetMesh
        Result from delaunay_tetrahedralization or constrained_delaunay_tetrahedralization
    scale : float
        Scale factor for each tetrahedron (0.0-1.0). Lower values create more separation.
    color_mode : str
        'random' - Random color per tetrahedron
        'quality' - Color by tetrahedron quality (red=bad, green=good)
        'index' - Color by tetrahedron index (gradient)
    min_volume : float
        Minimum volume threshold to skip degenerate tetrahedra
    show_edges : bool
        If True, show wireframe edges instead of solid faces

    Returns:
    --------
    trimesh.Trimesh or trimesh.Scene : The visualization object
    """
    tet_vertices = tet_mesh.vert_pos
    tet_indices = tet_mesh.tet_verts

    scaled_vertices = []
    vertex_colors = []
    num_tets = len(tet_indices)

    for i, tet in enumerate(tet_indices):
        # Get vertices for this tetrahedron
        v0 = np.array(tet_vertices[tet[0]])
        v1 = np.array(tet_vertices[tet[1]])
        v2 = np.array(tet_vertices[tet[2]])
        v3 = np.array(tet_vertices[tet[3]])

        # Calculate centroid
        centroid = (v0 + v1 + v2 + v3) / 4

        # Scale vertices around centroid
        v0_scaled = (v0 - centroid) * scale + centroid
        v1_scaled = (v1 - centroid) * scale + centroid
        v2_scaled = (v2 - centroid) * scale + centroid
        v3_scaled = (v3 - centroid) * scale + centroid

        # Check for degenerate tetrahedra
        edge1 = v1_scaled - v0_scaled
        edge2 = v2_scaled - v0_scaled
        edge3 = v3_scaled - v0_scaled
        volume = np.abs(np.dot(edge1, np.cross(edge2, edge3))) / 6.0

        if volume < min_volume:
            continue

        # Determine color based on mode
        if color_mode == 'random':
            color = np.random.rand(3)
        elif color_mode == 'quality':
            quality = compute_tet_quality(v0, v1, v2, v3)
            quality = np.clip(quality, 0, 1)
            # Red (bad) to Green (good)
            color = np.array([1.0 - quality, quality, 0.2])
        elif color_mode == 'index':
            t = i / max(num_tets - 1, 1)
            # Blue to Yellow gradient
            color = np.array([t, t, 1.0 - t])
        else:
            color = np.array([0.5, 0.5, 0.8])

        # Create four triangular faces for the tetrahedron
        # Face ordering for outward-facing normals
        scaled_vertices.extend([v0_scaled, v2_scaled, v1_scaled])  # Face 0
        scaled_vertices.extend([v0_scaled, v1_scaled, v3_scaled])  # Face 1
        scaled_vertices.extend([v1_scaled, v2_scaled, v3_scaled])  # Face 2
        scaled_vertices.extend([v0_scaled, v3_scaled, v2_scaled])  # Face 3

        # Add color for each vertex (12 vertices per tetrahedron)
        for _ in range(12):
            vertex_colors.append(color)

    if len(scaled_vertices) == 0:
        print("Warning: No valid tetrahedra to visualize")
        return None

    # Convert to numpy arrays
    scaled_vertices = np.array(scaled_vertices)
    vertex_colors = np.array(vertex_colors)

    # Create face indices
    face_indices = np.arange(len(scaled_vertices)).reshape(-1, 3)

    # Convert colors to RGBA (0-255)
    rgba_colors = np.zeros((len(vertex_colors), 4), dtype=np.uint8)
    rgba_colors[:, :3] = (vertex_colors * 255).astype(np.uint8)
    rgba_colors[:, 3] = 200  # Semi-transparent

    mesh = trimesh.Trimesh(
        vertices=scaled_vertices,
        faces=face_indices,
        vertex_colors=rgba_colors
    )

    if show_edges:
        # Create wireframe visualization
        edges = []
        for i in range(0, len(face_indices), 4):
            # Get the 4 faces of this tetrahedron
            base = i * 3
            # Extract unique edges from the 4 triangular faces
            tet_verts = scaled_vertices[base:base+12]
            if len(tet_verts) == 12:
                # 6 edges of a tetrahedron
                edges.extend([
                    [tet_verts[0], tet_verts[1]],
                    [tet_verts[0], tet_verts[2]],
                    [tet_verts[1], tet_verts[2]],
                    [tet_verts[0], tet_verts[5]],
                    [tet_verts[1], tet_verts[5]],
                    [tet_verts[2], tet_verts[5]],
                ])

        path = trimesh.load_path(edges)
        scene = trimesh.Scene([mesh, path])
        return scene

    return mesh


def visualize_with_surface(tet_mesh, manifold, tet_scale=0.7, surface_alpha=0.3):
    """
    Visualize tetrahedra alongside the original surface mesh.

    Parameters:
    -----------
    tet_mesh : TetMesh
        Result from constrained_delaunay_tetrahedralization
    manifold : Manifold
        The original manifold that was tetrahedralized
    tet_scale : float
        Scale factor for tetrahedra visualization
    surface_alpha : float
        Transparency of surface mesh (0-1)

    Returns:
    --------
    trimesh.Scene : Scene containing both visualizations
    """
    # Create tetrahedra visualization
    tet_vis = visualize_tetrahedra(tet_mesh, scale=tet_scale, color_mode='quality')

    # Create surface mesh visualization
    mesh = manifold.to_mesh()
    vertices = mesh.vert_properties[:, :3]
    faces = mesh.tri_verts

    surface_colors = np.zeros((len(vertices), 4), dtype=np.uint8)
    surface_colors[:, :3] = [100, 100, 200]  # Light blue
    surface_colors[:, 3] = int(surface_alpha * 255)

    surface_mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=faces,
        vertex_colors=surface_colors
    )

    # Combine into scene
    scene = trimesh.Scene()
    if tet_vis is not None:
        scene.add_geometry(tet_vis, geom_name='tetrahedra')
    scene.add_geometry(surface_mesh, geom_name='surface')

    return scene


def visualize_quality_histogram(tet_mesh):
    """
    Display a histogram of tetrahedron quality values.

    Parameters:
    -----------
    tet_mesh : TetMesh
        Result from tetrahedralization
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib required for histogram visualization")
        return

    tet_vertices = tet_mesh.vert_pos
    tet_indices = tet_mesh.tet_verts

    qualities = []
    for tet in tet_indices:
        v0 = np.array(tet_vertices[tet[0]])
        v1 = np.array(tet_vertices[tet[1]])
        v2 = np.array(tet_vertices[tet[2]])
        v3 = np.array(tet_vertices[tet[3]])
        q = compute_tet_quality(v0, v1, v2, v3)
        qualities.append(q)

    qualities = np.array(qualities)

    plt.figure(figsize=(10, 6))
    plt.hist(qualities, bins=50, edgecolor='black', alpha=0.7)
    plt.xlabel('Tetrahedron Quality')
    plt.ylabel('Count')
    plt.title(f'Tetrahedron Quality Distribution\n'
              f'Mean: {np.mean(qualities):.3f}, '
              f'Min: {np.min(qualities):.3f}, '
              f'Max: {np.max(qualities):.3f}')
    plt.axvline(np.mean(qualities), color='r', linestyle='--', label=f'Mean: {np.mean(qualities):.3f}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


def export_tetrahedra_to_obj(tet_mesh, filename, scale=1.0):
    """
    Export tetrahedra visualization to OBJ file.

    Parameters:
    -----------
    tet_mesh : TetMesh
        Result from tetrahedralization
    filename : str
        Output filename (should end in .obj)
    scale : float
        Scale factor for tetrahedra
    """
    mesh = visualize_tetrahedra(tet_mesh, scale=scale, color_mode='index')
    if mesh is not None:
        mesh.export(filename)
        print(f"Exported tetrahedra to {filename}")


def demo_visualization():
    """Demonstrate various visualization options."""
    print("=" * 60)
    print("TETRAHEDRALIZATION VISUALIZATION DEMO")
    print("=" * 60)

    # Create test shapes
    print("\n1. Creating test manifolds...")

    cube = m3d.Manifold.cube((2, 2, 2), center=True)
    sphere = m3d.Manifold.sphere(1.0, 24)

    # Boolean difference for more interesting shape
    shape = cube - sphere.scale((0.7, 0.7, 0.7))

    print(f"   Shape: {shape.num_vert()} vertices, {shape.num_tri()} triangles")

    # Tetrahedralize
    print("\n2. Performing constrained tetrahedralization...")
    tet_mesh = m3d.constrained_delaunay_tetrahedralization(shape, 0.0, 50)
    print(f"   Result: {tet_mesh.num_vert} vertices, {tet_mesh.num_tet} tetrahedra")

    # Show quality histogram
    print("\n3. Quality distribution:")
    tet_vertices = tet_mesh.vert_pos
    tet_indices = tet_mesh.tet_verts
    qualities = []
    for tet in tet_indices:
        v0 = np.array(tet_vertices[tet[0]])
        v1 = np.array(tet_vertices[tet[1]])
        v2 = np.array(tet_vertices[tet[2]])
        v3 = np.array(tet_vertices[tet[3]])
        qualities.append(compute_tet_quality(v0, v1, v2, v3))
    qualities = np.array(qualities)
    print(f"   Min quality: {np.min(qualities):.4f}")
    print(f"   Max quality: {np.max(qualities):.4f}")
    print(f"   Mean quality: {np.mean(qualities):.4f}")

    # Visualize
    print("\n4. Opening visualization...")
    print("   - Tetrahedra colored by quality (red=poor, green=good)")
    print("   - Surface mesh shown in blue")
    print("   - Close the window to continue")

    scene = visualize_with_surface(tet_mesh, shape, tet_scale=0.75, surface_alpha=0.2)
    scene.show()

    # Additional demos
    print("\n5. Additional visualization modes:")

    print("   a) Random colors (exploded view)...")
    vis = visualize_tetrahedra(tet_mesh, scale=0.6, color_mode='random')
    if vis:
        vis.show()

    print("   b) Index-based coloring...")
    vis = visualize_tetrahedra(tet_mesh, scale=0.8, color_mode='index')
    if vis:
        vis.show()

    # Simple example with random points
    print("\n6. Unconstrained tetrahedralization of random points...")
    np.random.seed(42)
    points = [(float(x), float(y), float(z)) for x, y, z in np.random.randn(50, 3)]
    tet_mesh_points = m3d.delaunay_tetrahedralization(points)
    print(f"   Result: {tet_mesh_points.num_tet} tetrahedra from {len(points)} points")

    vis = visualize_tetrahedra(tet_mesh_points, scale=0.7, color_mode='quality')
    if vis:
        vis.show()

    print("\nDemo complete!")


if __name__ == "__main__":
    demo_visualization()
