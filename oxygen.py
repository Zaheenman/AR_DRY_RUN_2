import numpy as np


C0 = 0.5
HD = 0.4
S_INITIAL = 0.97


def oxygen_consumption(M0, diameter, tissue_distance):

    r = diameter / 2
    rt = r + tissue_distance

    return np.pi * M0 * (rt**2 - r**2)


def saturation(x, x0, S0, Q, q):

    return S0 - (
        q / (Q * C0 * HD)
    ) * (x - x0)


def saturation_profile(
    S0,
    Q,
    diameter,
    length,
    M0,
    tissue_distance,
    points=100
):

    x = np.linspace(
        0,
        length,
        points
    )

    q = oxygen_consumption(
        M0,
        diameter,
        tissue_distance
    )

    S = saturation(
        x,
        0,
        S0,
        Q,
        q
    )

    return x, S
