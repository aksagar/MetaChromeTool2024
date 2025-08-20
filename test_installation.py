#!/usr/bin/env python3
"""
Test script to verify MetaChrome package installation.
"""

def test_imports():
    """Test that all modules can be imported."""
    try:
        print("Testing imports...")
        
        # Test main package import
        import metachrome
        print("✅ metachrome package imported successfully")
        
        # Test individual module imports
        from metachrome import ImageProcessor
        print("✅ ImageProcessor imported successfully")
        
        from metachrome import BatchProcessor
        print("✅ BatchProcessor imported successfully")
        
        from metachrome import SegmentationPostprocessing
        print("✅ SegmentationPostprocessing imported successfully")
        
        # Test version
        print(f"✅ MetaChrome version: {metachrome.__version__}")
        
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

def test_dependencies():
    """Test that required dependencies are available."""
    try:
        print("\nTesting dependencies...")
        
        import napari
        print("✅ napari imported successfully")
        
        import numpy
        print("✅ numpy imported successfully")
        
        import pandas
        print("✅ pandas imported successfully")
        
        import skimage
        print("✅ scikit-image imported successfully")
        
        import scipy
        print("✅ scipy imported successfully")
        
        import matplotlib
        print("✅ matplotlib imported successfully")
        
        import magicgui
        print("✅ magicgui imported successfully")
        
        import qtpy
        print("✅ qtpy imported successfully")
        
        import superqt
        print("✅ superqt imported successfully")
        
        import cellpose
        print("✅ cellpose imported successfully")
        
        import PIL
        print("✅ Pillow imported successfully")
        
        return True
        
    except ImportError as e:
        print(f"❌ Dependency import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Testing MetaChrome package installation...\n")
    
    # Test imports
    imports_ok = test_imports()
    
    # Test dependencies
    deps_ok = test_dependencies()
    
    print("\n" + "="*50)
    if imports_ok and deps_ok:
        print("🎉 All tests passed! MetaChrome is ready to use.")
        print("You can now run: metachrome")
    else:
        print("❌ Some tests failed. Please check the error messages above.")
        print("Make sure all dependencies are installed: pip install -r requirements.txt")

if __name__ == "__main__":
    main()
