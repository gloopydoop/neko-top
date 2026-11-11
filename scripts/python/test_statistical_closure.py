#!/usr/bin/env python3
"""Small dependency-light regression tests for the closure sign/conventions."""

import numpy as np

from statistical_closure import fit_boussinesq_nut, fit_scalar_gammat, reynolds_stress


def test_raw_to_centered():
    mean = np.array([[2.0], [-1.0], [0.5]])
    raw = np.zeros((6, 1))
    raw[:3, 0] = mean[:, 0] ** 2
    raw[3, 0] = mean[0, 0] * mean[1, 0]
    raw[4, 0] = mean[0, 0] * mean[2, 0]
    raw[5, 0] = mean[1, 0] * mean[2, 0]
    assert np.allclose(reynolds_stress(raw, mean), 0.0)


def test_manufactured_positive_nut():
    a = np.zeros((3, 3, 1))
    a[0, 1, 0] = 1.0
    a[1, 0, 0] = 1.0
    a[1, 2, 0] = -0.5
    a[2, 1, 0] = -0.5
    nut0 = 0.37
    r = -nut0 * a
    fitted, diagnostics = fit_boussinesq_nut(r, a, regularisation=1.0e-12, clip_negative=False)
    assert np.allclose(fitted, nut0, rtol=1.0e-10, atol=1.0e-12)
    assert diagnostics["mean_fit_relative_residual"] < 1.0e-10


def test_manufactured_positive_gammat():
    gradient = np.array([[2.0], [-1.0], [0.5]])
    gamma0 = 0.23
    flux = -gamma0 * gradient
    fitted, diagnostics = fit_scalar_gammat(
        flux, gradient, regularisation=1.0e-12, clip_negative=False
    )
    assert np.allclose(fitted, gamma0, rtol=1.0e-10, atol=1.0e-12)
    assert np.isclose(diagnostics["max_gammat"], gamma0)


if __name__ == "__main__":
    test_raw_to_centered()
    test_manufactured_positive_nut()
    test_manufactured_positive_gammat()
    print("statistical closure tests passed")
