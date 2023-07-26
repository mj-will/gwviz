#!/usr/bin/env python
"""Plot a diagram that shows the angles in the detector frame."""

from gwviz import plot_detector_angles
import numpy as np

# Define the parameters
phi = 0.7 * np.pi
theta = 0.2 * np.pi
psi = 0.4 * np.pi

# Units are in radians, so set degrees=False
fig = plot_detector_angles(phi, theta, psi)

fig.savefig("detector_angles.png")
