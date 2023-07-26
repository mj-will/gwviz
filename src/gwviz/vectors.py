"""Functions for creating vectors"""
import numpy as np


def get_basis_vectors(dims: int) -> np.ndarray:
    return np.eye(dims)


def angle_between_vectors(a: np.ndarray, b: np.ndarray) -> float:
    """Compute the angle in radians between two vectors."""
    return np.arccos(np.dot(a, b) / np.linalg.norm(a) * np.linalg.norm(b))
