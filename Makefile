.PHONY: help clean docs serve-docs install install-dev

# Default target
help:
	@echo "Available commands:"
	@echo "  install     - Install the package in development mode"
	@echo "  install-dev - Install with development dependencies"
	@echo "  docs        - Build the documentation"
	@echo "  serve-docs  - Build and serve documentation locally"
	@echo "  clean       - Clean build artifacts"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

docs:
	$(MAKE) -C docs html

serve-docs:
	$(MAKE) -C docs html
	cd docs/_build/html && python -m http.server 8000

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	$(MAKE) -C docs clean

# Documentation build targets
.PHONY: docs-clean docs-html

docs-clean:
	$(MAKE) -C docs clean

docs-html:
	$(MAKE) -C docs html
