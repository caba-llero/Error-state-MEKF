@echo off
echo Building and serving documentation...

REM Navigate to docs directory
cd docs

REM Build the documentation
echo Building documentation...
sphinx-build -b html . _build/html

REM Navigate to the built HTML directory
cd _build\html

REM Start the HTTP server
echo Starting server on port 8001...
echo Documentation will be available at: http://localhost:8001
echo Press Ctrl+C to stop the server
python -m http.server 8001
