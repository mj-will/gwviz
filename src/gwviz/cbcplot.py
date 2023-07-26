"""Functions for plotting CBC system with angles."""
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle
import numpy as np
from scipy.spatial.transform import Rotation as R

from .angles import to_radians
from .mplutils import add_arrow3d
from .orbit import (
    get_masses_coordinates,
    get_orbit_cartesian_coordinates,
    get_spins_coordinates,
    rotate_orbit,
)
from .vectors import angle_between_vectors, get_basis_vectors


def plot_cbc_diagram(
    s1,
    s2,
    L,
    iota,
    psi,
    phase,
    J=None,
    degrees=False,
    orbit_radius=1.0,
    n_orbit=100,
    axis_length=1.0,
    iota_arc_radius: float = 0.5,
    theta_jn_arc_radius: float = 0.5,
):
    """Plot a cbc diagram"""

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(vertical_axis="z", azim=45, elev=None)
    ax.set_axis_off()

    if degrees:
        iota = to_radians(iota)
        psi = to_radians(psi)
        phase = to_radians(phase)

    rotation = R.from_euler("yz", [iota, psi], degrees=False)

    s1 = np.array(s1)
    s2 = np.array(s2)
    L = np.array(L)

    if J is None:
        J = (s1 + s2) + L

    origin = np.zeros(3)
    bases = get_basis_vectors(3)

    ex_rot, ey_rot, ez_rot = rotation.apply(bases)

    arrow_kwargs = dict(
        shrinkA=0.0,
        shrinkB=0.0,
        arrowstyle=ArrowStyle.Fancy(head_length=2.0, head_width=2.0),
    )
    line_kwargs = dict(
        shrinkA=0.0,
        shrinkB=0.0,
        lw=0.5,
        arrowstyle="-",
        zorder=1,
    )

    angle_kwargs = dict(
        lw=1.0,
    )

    for b in bases:
        add_arrow3d(
            ax, origin, b, lw=1.0, zorder=0, color="grey", **arrow_kwargs
        )
    for b in [ex_rot, ey_rot]:
        add_arrow3d(
            ax, origin, b, lw=1.0, zorder=0, color="lightblue", **arrow_kwargs
        )

    # Orbits
    orbit = get_orbit_cartesian_coordinates(orbit_radius, n_orbit)
    rotated_orbit = rotate_orbit(rotation, orbit)

    ax.plot(*orbit, color="grey", lw=1.0)
    ax.plot(*rotated_orbit)

    # Masses
    m1_coords, m2_coords = get_masses_coordinates(
        phase,
        orbit_radius,
        rotation=rotation,
    )
    ax.scatter(*m1_coords, s=50.0, color="k")
    ax.scatter(*m2_coords, s=50.0, color="k")

    ax.text(*m1_coords, "$m_1$")
    ax.text(*m2_coords, "$m_2$")

    # Spins
    s1_coords, s2_coords = get_spins_coordinates(
        s1,
        s2,
        phase,
        orbit_radius,
        rotation=rotation,
    )
    add_arrow3d(ax, m1_coords, s1_coords, color="k", **arrow_kwargs)
    add_arrow3d(ax, m2_coords, s2_coords, color="k", **arrow_kwargs)
    ax.text(*s1_coords, r"$\mathbf{s_1}$")
    ax.text(*s2_coords, r"$\mathbf{s_2}$")

    # Momentum vectors
    L_rot = rotation.apply(L)
    J_rot = rotation.apply(J)

    add_arrow3d(ax, origin, L_rot, color="orange", **arrow_kwargs)
    J_rot = rotation.apply(J)
    add_arrow3d(ax, origin, J_rot, color="red", ls="-", **arrow_kwargs)
    ax.text(*J_rot, r"$\mathbf{J}$")
    ax.text(*L_rot, r"$\mathbf{L}$")

    # Angles
    L_angle_end = iota_arc_radius * L
    L_angle_end_trans = rotation.apply(L_angle_end)
    L_angle = angle_between_vectors(L_angle_end_trans, bases[2])

    J_angle_end = theta_jn_arc_radius * J
    J_angle_end_trans = rotation.apply(J_angle_end)

    theta_jn = angle_between_vectors(J_angle_end_trans, bases[2])

    J_rot = rotation.apply(J)
    J_angle = np.arctan2(J_rot[1], J_rot[0])

    rotation_J = R.from_euler("z", [J_angle], degrees=False)
    rotation_L = R.from_euler("z", [psi], degrees=False)

    theta = np.linspace(0, L_angle, 100)
    L_orbit = np.array(
        [
            iota_arc_radius * np.sin(theta),
            np.zeros(len(theta)),
            iota_arc_radius * np.cos(theta),
        ]
    )

    theta = np.linspace(0, theta_jn, 100)
    J_orbit = np.array(
        [
            theta_jn_arc_radius * np.sin(theta),
            np.zeros(len(theta)),
            theta_jn_arc_radius * np.cos(theta),
        ]
    )

    L_angle_line = rotation_L.apply(L_orbit.T)
    J_angle_line = rotation_J.apply(J_orbit.T)

    ax.plot(*L_angle_line.T, color="orange", **angle_kwargs)
    ax.plot(*J_angle_line.T, color="red", **angle_kwargs)

    ax.text(*L_angle_line[len(L_angle_line) // 2], r"$\iota$")
    ax.text(*J_angle_line[len(J_angle_line) // 2], r"$\theta_{JN}$")

    # Psi
    psi_line_rot = ex_rot.copy()
    psi_line_rot[2] = 0.0
    psi_line_rot /= np.linalg.norm(psi_line_rot)

    add_arrow3d(
        ax, origin, psi_line_rot, color="purple", ls=":", **line_kwargs
    )
    add_arrow3d(
        ax, psi_line_rot, ex_rot, color="purple", ls=":", **line_kwargs
    )

    theta = np.linspace(0, psi, 100)
    psi_r = 0.4
    psi_angle = np.array(
        [psi_r * np.cos(theta), psi_r * np.sin(theta), np.zeros(len(theta))]
    )
    ax.plot(*psi_angle, color="purple", **angle_kwargs)
    ax.text(*psi_angle.T[len(psi_angle.T) // 2], r"$\psi$")

    ax.text(*axis_length * bases[0], r"$\hat{\mathbf{x}}$")
    ax.text(*axis_length * bases[1], r"$\hat{\mathbf{y}}$")
    ax.text(*axis_length * bases[2], r"$\hat{\mathbf{n}}$")
    ax.text(*axis_length * ex_rot, r"$\hat{\mathbf{x}}'$")
    ax.text(*axis_length * ey_rot, r"$\hat{\mathbf{y}}'$")

    ax.set_xlim(axis_length * np.array([-1, 1]))
    ax.set_ylim(axis_length * np.array([-1, 1]))
    ax.set_zlim(axis_length * np.array([0.0, 1]))

    return fig
