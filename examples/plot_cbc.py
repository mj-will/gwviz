#!/usr/bin/env python
"""Plot a diagram of CBC system."""

from gwviz import plot_cbc_diagram
import numpy as np

# Define the parameters
s1 = [0.2, 0.3, 0.1]
s2 = [0.1, 0.0, 0.1]
L = [0, 0, 1.0]
phase = 1.3 * np.pi
iota = np.pi / 8
psi = 0.6 * np.pi

# Units are in radians, so set degrees=False
fig = plot_cbc_diagram(s1, s2, L, iota, psi, phase, degrees=False)

fig.savefig("cbc_diagram.png")
