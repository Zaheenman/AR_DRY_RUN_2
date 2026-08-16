import numpy as np


L0 = 0.01


def conducted_response(x, C, position):

    mask = x >= position

    y = x[mask]
    C_y = C[mask]

    weight = np.exp(
        -(y - position) / L0
    )

    return np.trapezoid(
        weight * C_y,
        y
    )


def conducted_profile(x, C):

    return np.array([
        conducted_response(
            x,
            C,
            position
        )
        for position in x
    ])
