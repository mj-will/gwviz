import matplotlib.pyplot as plt
import pytest


@pytest.fixture()
def close_figures():
    """Function to close all open figures"""

    def fn():
        plt.close("all")

    return fn
