"""Functions related to angles."""
import numpy as np


def to_radians(angle: float) -> float:
    """Convert an angle from degrees to radians."""
    return angle * np.pi / 180
