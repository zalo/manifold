#!/bin/bash
# Deploy tetrahedralization visualization to Cloudflare Pages
#
# Usage: ./deploy.sh
#
# Prerequisites:
#   - wrangler installed: npm install -g wrangler
#   - Cloudflare authentication: wrangler login
#   - Python with manifold3d built

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PROJECT_NAME="manifold-tet-viewer"
MAIN_URL="https://${PROJECT_NAME}.pages.dev"

echo "========================================"
echo "Manifold Tetrahedralization Deployment"
echo "========================================"

# Step 1: Run tests and export data
echo ""
echo "Step 1: Exporting test data..."
python3 export_test_data.py

if [ ! -f "public/test_data.json" ]; then
    echo "Error: test_data.json was not created"
    exit 1
fi

echo ""
echo "Step 2: Deploying to Cloudflare Pages..."

# Check if wrangler is installed
if ! command -v wrangler &> /dev/null; then
    echo "Error: wrangler is not installed"
    echo "Install with: npm install -g wrangler"
    exit 1
fi

# Create project if it doesn't exist
wrangler pages project create "$PROJECT_NAME" --production-branch main 2>/dev/null || true

# Deploy using wrangler pages
wrangler pages deploy public --project-name="$PROJECT_NAME" --commit-dirty=true

echo ""
echo "========================================"
echo "Deployment complete!"
echo ""
echo "Main URL: $MAIN_URL"
echo "========================================"
