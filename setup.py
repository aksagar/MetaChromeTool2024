from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="metachrome",
    version="1.0.0",
    author="Md Abdul Kader Sagar",
    author_email="sagarm2@nih.gov",
    description="Metaphase Chromosome Analysis with Napari",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/metachrome",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Image Processing",
    ],
    python_requires=">=3.8",
    install_requires=[
        "napari>=0.4.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scikit-image>=0.18.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "magicgui>=0.2.0",
        "qtpy>=1.9.0",
        "superqt>=0.2.0",
        "cellpose>=2.0.0",
        "Pillow>=8.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
    },
    entry_points={
        "console_scripts": [
            "metachrome=metachrome.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "metachrome": ["*.py"],
    },
)
