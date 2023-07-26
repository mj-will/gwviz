from gwviz.cbcplot import plot_cbc_diagram
import pytest


@pytest.fixture()
def parameters():
    return (
        [0.2, 0.3, 0.1],
        [0.1, 0.0, 0.1],
        [0, 0, 1.0],
        0.4,
        0.0,
        0.0,
    )


def test_defaults(parameters, close_figures):
    plot_cbc_diagram(*parameters)
    close_figures()
