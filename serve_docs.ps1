Write-Host "Building and serving documentation..." -ForegroundColor Green

# Navigate to docs directory
Set-Location docs

# Build the documentation
Write-Host "Building documentation..." -ForegroundColor Yellow
sphinx-build -b html . _build/html

# Navigate to the built HTML directory
Set-Location _build/html

# Start the HTTP server
Write-Host "Starting server on port 8001..." -ForegroundColor Yellow
Write-Host "Documentation will be available at: http://localhost:8001" -ForegroundColor Cyan
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Magenta
python -m http.server 8001
