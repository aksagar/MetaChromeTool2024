"""
Metaphase Chromosome Analysis with Napari

A Python package for analyzing metaphase chromosomes using the Napari platform.
It facilitates the visualization and segmentation of chromosome images, enabling users 
to efficiently assess chromosome structures and perform quantitative analysis.
The code integrates tools for detecting centromeres and measuring CENP-A levels 
within metaphase chromosome regions, enhancing the accuracy of chromosome analysis.

Author: Md Abdul Kader Sagar
Email: sagarm2@nih.gov
Institute: National Cancer Institute/NIH
"""

__version__ = "1.0.0"
__author__ = "Md Abdul Kader Sagar"
__email__ = "sagarm2@nih.gov"
__institute__ = "National Cancer Institute/NIH"

from .main import *
from .image_processor import ImageProcessor
from .batch_processor import BatchProcessor
from .segmentation_postprocessing import SegmentationPostprocessing

__all__ = [
    "ImageProcessor",
    "BatchProcessor", 
    "SegmentationPostprocessing",
]
