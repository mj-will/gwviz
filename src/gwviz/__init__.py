"""GWViz: Visualizing gravitational-wave systems using Python"""
from importlib.metadata import PackageNotFoundError, version

from .cbcplot import plot_cbc_diagram
from .detplot import plot_detector_angles

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"

__all__ = [
    "plot_cbc_diagram",
    "plot_detector_angles",
]
