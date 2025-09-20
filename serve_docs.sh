#!/bin/bash
echo "Building and serving documentation..."

# Navigate to docs directory
cd docs

# Build the documentation
echo "Building documentation..."
sphinx-build -b html . _build/html

# Navigate to the built HTML directory
cd _build/html

# Start the HTTP server
echo "Starting server on port 8001..."
echo "Documentation will be available at: http://localhost:8001"
echo "Press Ctrl+C to stop the server"
python -m http.server 8001
