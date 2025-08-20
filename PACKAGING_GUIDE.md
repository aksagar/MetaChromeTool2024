# MetaChrome Packaging Guide

This document explains how the MetaChrome project has been packaged for pip installation.

## Package Structure

```
MetaChromeTool2024/
├── metachrome/                    # Main package directory
│   ├── __init__.py               # Package initialization
│   ├── main.py                   # Main application entry point
│   ├── image_processor.py        # Image processing functionality
│   ├── batch_processor.py        # Batch processing functionality
│   └── segmentation_postprocessing.py  # Segmentation post-processing
├── setup.py                      # Traditional setup configuration
├── pyproject.toml                # Modern Python packaging configuration
├── requirements.txt              # Python dependencies
├── MANIFEST.in                   # Files to include in distribution
├── LICENSE                       # MIT License
├── install.py                    # Installation helper script
├── test_installation.py          # Installation verification script
└── README.md                     # Project documentation
```

## Installation Methods

### Method 1: pip install (when published to PyPI)
```bash
pip install metachrome
```

### Method 2: Local development installation
```bash
# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

### Method 3: Using the installation script
```bash
python install.py
```

## Usage

After installation, you can use MetaChrome in several ways:

### Command Line Interface
```bash
metachrome
```

### Python Import
```python
import metachrome
from metachrome import ImageProcessor, BatchProcessor, SegmentationPostprocessing

# Use the classes
processor = ImageProcessor()
```

### Direct Script Execution
```bash
python -m metachrome.main
```

## Testing Installation

Run the test script to verify everything is working:
```bash
python test_installation.py
```

## Package Configuration

### Entry Points
The package defines a console script entry point:
- `metachrome=metachrome.main:main`

### Dependencies
All required dependencies are specified in:
- `requirements.txt` - for manual installation
- `setup.py` and `pyproject.toml` - for pip installation

### Version Management
- Version is set to "1.0.0" in `__init__.py`
- Can be updated to use dynamic versioning with setuptools_scm

## Building for Distribution

To build the package for distribution:

```bash
# Install build tools
pip install build

# Build the package
python -m build

# This creates dist/ directory with wheel and source distribution
```

## Publishing to PyPI

1. Build the package: `python -m build`
2. Upload to PyPI: `python -m twine upload dist/*`
3. Users can then install with: `pip install metachrome`

## Development

For development, install in editable mode:
```bash
pip install -e .
```

This allows you to modify the code and see changes immediately without reinstalling.

## Notes

- The original code structure has been preserved
- All Python files are copied to the `metachrome/` package directory
- The main entry point is `metachrome.main:main`
- Dependencies are carefully specified to ensure compatibility
- The package follows Python packaging best practices
