import numpy as np
import pytest

from sat_trace.sat_trace import SatelliteTrace, plot_traces


def test_init_sma_normal():
    s = SatelliteTrace(40708, 0.8320, 61, 120, 270)
    assert s.sma == 40708


def test_init_sma_low():
    with pytest.raises(Exception):
        SatelliteTrace(100, 0.8320, 61, 120, 270)


def test_init_sma_high():
    with pytest.raises(Exception):
        SatelliteTrace(10000000, 0.8320, 61, 120, 270)


def test_init_ecc_high():
    with pytest.raises(Exception):
        SatelliteTrace(100, 1, 61, 120, 270)


def test_init_ecc_low():
    with pytest.raises(Exception):
        SatelliteTrace(100, -1, 61, 120, 270)


def test_init_radius_true():
    s = SatelliteTrace(40708, 0.8320, 61, 120, 270)
    assert round(s.radius_ap) == 74577


def test_correction_for_t():
    s = SatelliteTrace(40708, 0.8320, 61, 120, 270)
    correction, factor = s._determine_correction_for_t(
        np.radians(146), np.radians(-270)
    )
    assert round(np.degrees(correction)) == round(np.degrees(-2 * np.pi))
    assert factor == 1


def test_correction_for_longitude():
    s = SatelliteTrace(40708, 0.8320, 61, 120, 270)
    correction, factor = s._determine_correction_for_lo(
        np.radians(270), np.radians(30)
    )
    assert round(np.degrees(correction)) == round(np.degrees(2 * np.pi))
    assert factor == +1


# def test_plot_traces_good_input():
#     s1 = SatelliteTrace(40708, 0.8320, 61, 120, 270)
#     s2 = SatelliteTrace(40708, 0.8320, 61, 120, 270)
#     s3 = SatelliteTrace(40708, 0.8320, 61, 120, 270)
#     plot_traces(s1, s2, s3)


def test_plot_traces_no_input():
    with pytest.raises(Exception):
        plot_traces()
