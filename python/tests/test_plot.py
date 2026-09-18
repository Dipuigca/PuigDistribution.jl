"""Smoke-test de la visualización (matplotlib backend Agg)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from puigdist.core import PuigDistribution, MB
from puigdist.plot import Puig_plot


def test_puig_plot_all():
    d = PuigDistribution(3.0, 2.0, 1.0)
    fig = Puig_plot(d, options="all")
    assert fig is not None
    assert len(fig.axes) == 4
    plt.close(fig)


def test_puig_plot_tupla():
    d = PuigDistribution(3.0, 2.0, 1.0)
    fig = Puig_plot(d, options=("pdf", "surv"))
    assert len(fig.axes) == 2
    plt.close(fig)


def test_puig_plot_single():
    d = PuigDistribution(3.0, 2.0, 1.0)
    for o in ("pdf", "cdf", "surv", "data"):
        fig = Puig_plot(d, options=o)
        assert fig is not None
        plt.close(fig)


def test_puig_plot_mb():
    mb = MB(2.0, 1.0)
    fig = Puig_plot(mb, options=("pdf", "surv"))
    assert len(fig.axes) == 2
    plt.close(fig)


def test_puig_plot_bad_option():
    import pytest

    d = PuigDistribution(3.0, 2.0, 1.0)
    with pytest.raises(ValueError):
        Puig_plot(d, options="xxx")