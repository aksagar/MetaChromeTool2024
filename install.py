#!/usr/bin/env python3
"""
Installation script for MetaChrome package.
This script installs the MetaChrome package in development mode.
"""

import subprocess
import sys
import os

def install_package():
    """Install the MetaChrome package in development mode."""
    try:
        # Install in development mode
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "."])
        print("✅ MetaChrome package installed successfully!")
        print("You can now run: metachrome")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error installing package: {e}")
        return False

def install_dependencies():
    """Install required dependencies."""
    try:
        print("Installing dependencies...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
        print("✅ Dependencies installed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error installing dependencies: {e}")
        return False

def main():
    """Main installation function."""
    print("🚀 Installing MetaChrome package...")
    
    # Install dependencies first
    if not install_dependencies():
        print("Failed to install dependencies. Please check your Python environment.")
        return
    
    # Install the package
    if install_package():
        print("\n🎉 Installation complete!")
        print("You can now use MetaChrome by running:")
        print("  metachrome")
        print("\nOr import it in Python:")
        print("  import metachrome")
    else:
        print("\n❌ Installation failed. Please check the error messages above.")

if __name__ == "__main__":
    main()
