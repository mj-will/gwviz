"""Functions for plotting the detector angles"""
from typing import Dict, List, Optional

from matplotlib.figure import Figure
from matplotlib.patches import ArrowStyle
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R

from .angles import to_radians
from .mplutils import add_arrow3d
from .vectors import get_basis_vectors


def plot_detector_angles(
    phi: float = 2.2,
    theta: float = 0.63,
    psi: float = 1.26,
    degrees: bool = False,
    n_angle: int = 100,
    distance: float = 2.0,
    source_frame_labels: Optional[List[str]] = None,
    source_frame_labels_sky: Optional[List[str]] = None,
    detector_frame_labels: Optional[List[str]] = None,
    xlims: Optional[List[float]] = None,
    ylims: Optional[List[float]] = None,
    zlims: Optional[List[float]] = None,
    view_kwargs: Optional[Dict] = None,
    phi_label_position: float = 0.5,
    theta_label_position: float = 0.5,
    psi_label_position: float = 0.5,
    phi_arc_radius: float = 0.5,
    theta_arc_radius: float = 0.5,
    psi_arc_radius: float = 0.5,
    phi_label_pad: Optional[List[float]] = None,
    theta_label_pad: Optional[List[float]] = None,
    psi_label_pad: Optional[List[float]] = None,
    label_pad: Optional[List[float]] = None,
) -> Figure:
    """Plot the detector angles.

    Parameters
    ----------
    phi :
        Polar angle about the z-axis.
    theta :
        Polar angle from the z-axis.
    psi :
        Polarization angle.
    degrees :
        If True, angles are assumed to be in degrees. If False, angles are
        assumed to be in radians.
    n_angle :
        Number of samples to use when plotting the angles (arcs).
    view_kwargs :
        Keyword arguments passed to `view_init`.
    phi_label_position :
        Relative position of the label, 0 corresponds to the start of the arc
        and 1 to the end.
    theta_label_position :
        Relative position of the label, 0 corresponds to the start of the arc
        and 1 to the end.
    psi_label_position :
        Relative position of the label, 0 corresponds to the start of the arc
        and 1 to the end.
    phi_arc_radius :
        Radius of the arc plot to show the phi angle.
    theta_arc_radius :
        Radius of the arc plot to show the theta angle.
    psi_arc_radius :
        Radius of the arc plot to show the psi angle.
    """
    fig = plt.figure()

    if degrees:
        phi = to_radians(phi)
        theta = to_radians(theta)
        psi = to_radians(psi)

    if phi_label_pad is None:
        phi_label_pad = [0.0, 0.0, -0.1]
    if theta_label_pad is None:
        theta_label_pad = [0.0, 0.0, 0.1]
    if psi_label_pad is None:
        psi_label_pad = [0.0, 0.0, 0.1]
    if label_pad is None:
        label_pad = [0.0, 0.0, 0.1]

    phi_label_pad = np.array(phi_label_pad)
    theta_label_pad = np.array(theta_label_pad)
    psi_label_pad = np.array(psi_label_pad)
    label_pad = np.array(label_pad)

    # Define the rotation matrices
    rotation_sky = R.from_euler("yz", [theta, phi], degrees=False)
    rotation_psi = R.from_euler("z", [-psi], degrees=False)

    fig = plt.figure()

    ax = fig.add_subplot(111, projection="3d")

    default_view_kwargs = dict(
        vertical_axis="z",
        azim=45,
        elev=20,
    )
    if view_kwargs is not None:
        default_view_kwargs.update(view_kwargs)

    ax.view_init(**default_view_kwargs)
    ax.set_axis_off()

    origin = np.zeros(3)
    detector_frame_axes = get_basis_vectors(3)

    # Rotate the vector that points along the LOS
    source_frame_origin_unit = rotation_sky.apply(detector_frame_axes[-1])
    source_frame_origin = distance * source_frame_origin_unit

    # Define the source frame axes
    source_frame_axes = detector_frame_axes.copy()
    source_frame_axes[[0, 2]] *= -1

    # Rotate the source frame
    # With psi
    rotated_source_frame_axes = rotation_sky.apply(
        rotation_psi.apply(source_frame_axes)
    )
    # Without psi
    source_frame_axes_sky = rotation_sky.apply(source_frame_axes)

    # Configure the axes arrows
    axes_arrows_kwargs = dict(
        arrowstyle=ArrowStyle.Fancy(head_length=2.0, head_width=2.0),
    )

    # Labels for the different axes
    if detector_frame_labels is None:
        detector_frame_labels = (r"$\hat{p}$", r"$\hat{q}$", r"$\hat{z}$")
    if source_frame_labels is None:
        source_frame_labels = (r"$e_1$", r"$e_2$", r"$\hat{n}$")
    if source_frame_labels_sky is None:
        source_frame_labels_sky = (r"$e_1^s$", r"$e_2^s$")

    # Plot the detector frame axes
    for vec, label in zip(detector_frame_axes, detector_frame_labels):
        add_arrow3d(ax, origin, vec, color="k", **axes_arrows_kwargs)
        ax.text(*vec + label_pad, label, ha="center", va="center")

    # Add the arrow pointing to the source
    add_arrow3d(ax, origin, source_frame_origin, color="grey", ls=":")

    # Plot the rotated source frame
    for vec, label in zip(rotated_source_frame_axes, source_frame_labels):
        vec = source_frame_origin + vec
        add_arrow3d(
            ax, source_frame_origin, vec, color="k", **axes_arrows_kwargs
        )
        ax.text(*vec + label_pad, label, ha="center", va="center")

    # Plot the source frame without the psi rotation
    for vec, label in zip(source_frame_axes_sky[:2], source_frame_labels_sky):
        vec = source_frame_origin + vec
        add_arrow3d(
            ax,
            source_frame_origin,
            vec,
            color="grey",
            **axes_arrows_kwargs,
        )
        ax.text(*vec + label_pad, label)

    # Plot the projections for the line of sight
    los_projection_xy = source_frame_origin.copy()
    los_projection_xy[2] = 0.0

    add_arrow3d(ax, origin, los_projection_xy, color="grey", ls=":")
    add_arrow3d(
        ax, source_frame_origin, los_projection_xy, color="grey", ls=":"
    )

    # Plot the angle arcs
    phi_vec = np.linspace(0, phi, n_angle)
    phi_angle = np.array(
        [
            phi_arc_radius * np.cos(phi_vec),
            phi_arc_radius * np.sin(phi_vec),
            np.zeros(len(phi_vec)),
        ]
    )
    ax.plot(*phi_angle, color="blue")

    ax.text(
        *phi_angle[:, int(phi_label_position * n_angle)] + phi_label_pad,
        r"$\phi$",
        verticalalignment="center",
        horizontalalignment="center",
    )

    theta_vec = np.linspace(0, theta, n_angle)
    theta_angle = np.array(
        [
            theta_arc_radius * np.sin(theta_vec),
            np.zeros(len(theta_vec)),
            theta_arc_radius * np.cos(-theta_vec),
        ]
    )

    theta_angle_rot = (
        R.from_euler(
            "z",
            [
                phi,
            ],
        )
        .apply(theta_angle.T)
        .T
    )

    ax.plot(*theta_angle_rot, color="turquoise")
    ax.text(
        *theta_angle_rot[:, int(theta_label_position * n_angle)]
        + theta_label_pad,
        r"$\theta$",
        verticalalignment="center",
        horizontalalignment="center",
    )

    psi_vec = np.linspace(0, -psi, n_angle) + np.pi
    psi_angle = np.array(
        [
            psi_arc_radius * np.cos(psi_vec),
            psi_arc_radius * np.sin(psi_vec),
            np.zeros(len(psi_vec)),
        ]
    )
    psi_angle_rot = (
        source_frame_origin[:, np.newaxis] + rotation_sky.apply(psi_angle.T).T
    )
    ax.plot(*psi_angle_rot, color="purple")
    ax.text(
        *psi_angle_rot[:, int(psi_label_position * n_angle)] + psi_label_pad,
        r"$\psi$",
        verticalalignment="center",
        horizontalalignment="center",
    )

    if xlims is None:
        xlims = [-1, 1]
    if ylims is None:
        ylims = [0, 1.3]
    if zlims is None:
        zlims = [0, 2]

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_zlim(zlims)

    ax.set_box_aspect((np.ptp(xlims), np.ptp(ylims), np.ptp(zlims)), zoom=1.0)

    return fig
