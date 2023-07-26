"""Functions related to plotting orbits"""
import numpy as np


def get_orbit_cartesian_coordinates(radius, n):
    """Get the Cartesian coordinates of an orbit in the x-y plane."""
    theta = np.linspace(0, 2 * np.pi, n)
    orbit = np.array(
        [
            radius * np.cos(theta),
            radius * np.sin(theta),
            np.zeros(len(theta)),
        ]
    )
    return orbit


def rotate_orbit(rotation, orbit):
    """Rotate an orbit"""
    return rotation.apply(orbit.T).T


def get_masses_coordinates(angle, radius, rotation=None):
    """Get the coordinates for the masses"""
    m1_coords = np.array([radius * np.cos(angle), radius * np.sin(angle), 0.0])
    m2_coords = np.array(
        [radius * np.cos(angle + np.pi), radius * np.sin(angle + np.pi), 0.0]
    )

    if rotation is not None:
        m1_coords = rotation.apply(m1_coords)
        m2_coords = rotation.apply(m2_coords)

    return m1_coords, m2_coords


def get_spins_coordinates(s1, s2, angle, radius, rotation=None):
    """Get the coordinates of the end of the spin arrows"""
    m1_coords, m2_coords = get_masses_coordinates(angle, radius)
    s1_coords = m1_coords + s1
    s2_coords = m2_coords + s2

    if rotation is not None:
        s1_coords = rotation.apply(s1_coords)
        s2_coords = rotation.apply(s2_coords)

    return s1_coords, s2_coords
