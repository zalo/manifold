# Manifold Tetrahedralization Viewer

Interactive Three.js visualization for Delaunay tetrahedralization test results.

**Live Demo:** https://manifold-tet-viewer.pages.dev

## Features

- View tetrahedralization results from multiple test cases
- Color-coded quality visualization (red = poor, green = good)
- Adjustable tetrahedron scale (exploded view)
- Optional surface mesh overlay
- Wireframe mode
- Auto-rotation
- Statistics panel (vertices, tetrahedra, quality metrics)

## Prerequisites

- Python 3.9+ with manifold3d built (with tetrahedralization support)
- Node.js with wrangler: `npm install -g wrangler`
- Cloudflare account and authentication: `wrangler login`

## Usage

### Local Development

```bash
# Generate test data
python3 export_test_data.py

# Start local server
python3 -m http.server 8080 -d public

# Open http://localhost:8080 in browser
```

### Deploy to Cloudflare

```bash
# One-command deploy (exports data + deploys)
./deploy.sh
```

Or using npm:

```bash
npm run deploy
```

## Test Cases

The visualization includes 10 test cases:

1. **Random Points** - Unconstrained Delaunay of 50 random points
2. **Cube** - Constrained tetrahedralization of a cube
3. **Sphere** - Constrained tetrahedralization of a geodesic sphere
4. **Cylinder** - Constrained tetrahedralization of a cylinder
5. **Torus** - Constrained tetrahedralization of a torus
6. **Boolean Difference** - Cube minus sphere
7. **Extruded L-Shape** - Extruded L-shaped cross section
8. **High-Genus Shape** - Cube with cylindrical holes (genus 3)
9. **Tetrahedron** - Single tetrahedron (perfect quality)
10. **Large Point Cloud** - Unconstrained Delaunay of 200 points

## Quality Metric

Tetrahedron quality is computed as:

```
quality = (12 / sqrt(2)) * volume / rms_edge_length^3
```

Where:
- `volume` is the signed tetrahedron volume
- `rms_edge_length` is the root mean square of all 6 edge lengths
- A regular tetrahedron has quality = 1.0
- Degenerate (flat) tetrahedra have quality = 0.0

## Files

- `export_test_data.py` - Generates test data JSON from manifold3d
- `public/index.html` - Three.js visualization webpage
- `public/test_data.json` - Generated test data (not committed)
- `wrangler.toml` - Cloudflare Pages configuration
- `deploy.sh` - Deployment script
