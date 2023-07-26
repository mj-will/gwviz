from gwviz.detplot import plot_detector_angles


def test_defaults(close_figures):
    """Test plotting with the defaults"""
    plot_detector_angles()
    close_figures()
