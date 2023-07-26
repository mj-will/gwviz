"""Utilities for matplotlib"""
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np


class Arrow3D(FancyArrowPatch):
    """Arrow class the works in 3D plots.

    Based on: https://stackoverflow.com/a/74122407
    """

    def __init__(self, start, end, shrinkA=0, shrinkB=0, *args, **kwargs):
        xs = np.array([start[0], end[0]])
        ys = np.array([start[1], end[1]])
        zs = np.array([start[2], end[2]])
        FancyArrowPatch.__init__(
            self,
            (0, 0),
            (0, 0),
            shrinkA=shrinkA,
            shrinkB=shrinkB,
            *args,
            **kwargs,
        )
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def add_arrow3d(ax, start, end, *args, **kwargs):
    """Add an arrow to a 3d plot"""
    arrow = Arrow3D(start, end, *args, **kwargs)
    ax.add_artist(arrow)
