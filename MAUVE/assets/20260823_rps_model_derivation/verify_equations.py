#!/usr/bin/env python3
"""Numerically falsify key algebraic identities in the RPS derivation.

This uses only the Python standard library so the verification is reproducible
in the current workspace without installing a scientific runtime.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


OUT = Path(__file__).resolve().parent / "verification_results.json"


def central_difference(func, x: float, step: float) -> float:
    return (func(x + step) - func(x - step)) / (2.0 * step)


def relerr(actual: float, expected: float) -> float:
    scale = max(abs(expected), 1e-15)
    return abs(actual - expected) / scale


def adaptive_simpson(func, a: float, b: float, eps: float = 1e-11, depth: int = 22) -> float:
    def simp(x0: float, x1: float) -> float:
        xm = 0.5 * (x0 + x1)
        return (x1 - x0) * (func(x0) + 4.0 * func(xm) + func(x1)) / 6.0

    whole = simp(a, b)

    def rec(x0: float, x1: float, value: float, tol: float, left: int) -> float:
        xm = 0.5 * (x0 + x1)
        vl, vr = simp(x0, xm), simp(xm, x1)
        if left == 0 or abs(vl + vr - value) <= 15.0 * tol:
            return vl + vr + (vl + vr - value) / 15.0
        return rec(x0, xm, vl, tol / 2.0, left - 1) + rec(xm, x1, vr, tol / 2.0, left - 1)

    return rec(a, b, whole, eps, depth)


def verify_miyamoto_nagai() -> dict:
    g, mass, a, b = 1.7, 3.2, 1.3, 0.27
    errors_r, errors_z = [], []

    def phi(r: float, z: float) -> float:
        big_b = math.sqrt(z * z + b * b)
        return -g * mass / math.sqrt(r * r + (a + big_b) ** 2)

    for r, z in [(0.8, 0.11), (2.1, 0.7), (4.2, -0.45), (7.0, 1.4)]:
        big_b = math.sqrt(z * z + b * b)
        d = math.sqrt(r * r + (a + big_b) ** 2)
        analytic_r = g * mass * r / d**3
        analytic_z = g * mass * (a + big_b) * z / (big_b * d**3)
        numeric_r = central_difference(lambda rr: phi(rr, z), r, 1e-5 * max(1.0, r))
        numeric_z = central_difference(lambda zz: phi(r, zz), z, 1e-5 * max(1.0, abs(z)))
        errors_r.append(relerr(numeric_r, analytic_r))
        errors_z.append(relerr(numeric_z, analytic_z))
    return {"max_relative_error_dPhi_dR": max(errors_r), "max_relative_error_dPhi_dz": max(errors_z)}


def verify_nfw() -> dict:
    g, m200, rs, c = 2.3, 8.1, 1.7, 5.2
    ac = math.log1p(c) - c / (1.0 + c)

    def phi(r: float) -> float:
        return -g * m200 * math.log1p(r / rs) / (r * ac)

    errors = []
    for r in [0.13, 0.7, 2.4, 8.0, 21.0]:
        x = r / rs
        enclosed = m200 * (math.log1p(x) - x / (1.0 + x)) / ac
        expected = g * enclosed / r**2
        numeric = central_difference(phi, r, 1e-6 * max(1.0, r))
        errors.append(relerr(numeric, expected))
    return {"max_relative_error_force": max(errors)}


def verify_edge_on_integral() -> dict:
    # Compare direct x integration with the independently transformed cosh integral.
    values = []
    for b_over_rg in [0.05, 0.2, 0.7, 1.5, 3.5]:
        b = b_over_rg
        direct = 2.0 * adaptive_simpson(
            lambda x: math.exp(-math.sqrt(x * x + b * b)), 0.0, 40.0, eps=2e-10
        )
        transformed = 2.0 * b * adaptive_simpson(
            lambda t: math.exp(-b * math.cosh(t)) * math.cosh(t), 0.0, 14.0, eps=2e-10
        )
        values.append(relerr(direct, transformed))
    # The full direct central chord is integral exp(-|x|) dx = 2.
    central = 2.0 * adaptive_simpson(lambda x: math.exp(-x), 0.0, 40.0)
    return {
        "max_relative_error_direct_vs_bessel_representation": max(values),
        "central_chord": central,
        "central_chord_expected": 2.0,
        "central_relative_error": relerr(central, 2.0),
    }


def verify_point_mass_extrema() -> dict:
    # Dimensionless R0=GM=1.
    n = 400_001
    zmax, gzmax = 0.0, -1.0
    for i in range(n):
        z = 3.0 * i / (n - 1)
        gz = z / (1.0 + z * z) ** 1.5
        if gz > gzmax:
            zmax, gzmax = z, gz
    n2 = 400_001
    rmax, grmax = 1.0, -1.0
    for i in range(n2):
        r = 1.0 + 4.0 * i / (n2 - 1)
        gr = 1.0 / r**2 - 1.0 / r**3
        if gr > grmax:
            rmax, grmax = r, gr
    return {
        "vertical_z_at_max": zmax,
        "vertical_z_expected": 1.0 / math.sqrt(2.0),
        "vertical_gmax": gzmax,
        "vertical_gmax_expected": 2.0 / (3.0 * math.sqrt(3.0)),
        "radial_R_at_max": rmax,
        "radial_R_expected": 1.5,
        "radial_gmax": grmax,
        "radial_gmax_expected": 4.0 / 27.0,
    }


def verify_impulse_limits() -> dict:
    vc, vesc = 200.0, 500.0

    def jcrit(wphi: float) -> float:
        return -vc * wphi + math.sqrt(vc * vc * wphi * wphi + vesc * vesc - vc * vc)

    return {
        "wphi_0": jcrit(0.0),
        "wphi_0_expected": math.sqrt(vesc * vesc - vc * vc),
        "wphi_plus1": jcrit(1.0),
        "wphi_plus1_expected": vesc - vc,
        "wphi_minus1": jcrit(-1.0),
        "wphi_minus1_expected": vesc + vc,
    }


def main() -> None:
    results = {
        "miyamoto_nagai": verify_miyamoto_nagai(),
        "nfw": verify_nfw(),
        "edge_on_integral": verify_edge_on_integral(),
        "point_mass_extrema": verify_point_mass_extrema(),
        "impulse_limits_km_s": verify_impulse_limits(),
    }
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
