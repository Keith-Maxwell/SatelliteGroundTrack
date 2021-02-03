import numpy as np
from sat_trace.sat_trace import Satellite

s = Satellite(40708, 0.8320, 61, 120, 270)


def test_init():
    assert round(s.radius_ap) == 74577


def test_correction_for_t():
    correction, factor = s._determine_correction_for_t(np.radians(146), np.radians(-270))
    assert round(np.degrees(correction)) == round(np.degrees(-2 * np.pi))
    assert factor == 1


def test_correction_for_longitude():
    # throws error, don't know why. It works on the real case, which is the same
    correction, factor = s._determine_correction_for_lo(np.radians(270), np.radians(30))
    assert round(np.degrees(correction)) == round(np.degrees(np.pi))
    assert factor == -1


if __name__ == "__main__":
    test_init()
    test_correction_for_t()
    test_correction_for_longitude()
    print("ok")
