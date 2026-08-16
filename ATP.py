import numpy as np


HD = 0.4
HT = 0.3
C0_OXYGEN = 0.5

R0 = 1.4e-3
R1 = 0.891
KD = 2e-6

C_INITIAL = 0.5e-3


def atp_release(S):

    return R0 * (1 - R1 * S)


def alpha(D, S0, q):

    return (
        HT * R0 / (4 * KD)
        * (
            D * (1 - R1 * S0)
            - (
                (1 - HD) * R1 * q
                / (
                    np.pi
                    * C0_OXYGEN
                    * HD
                    * KD
                )
            )
        )
    )


def beta(D, Q, q):

    return (
        D * HT * R0 * R1 * q
        / (
            4
            * Q
            * C0_OXYGEN
            * HD
            * KD
        )
    )


def gamma(D, Q):

    if not np.isfinite(D) or D <= 0:
        raise ValueError(
            f"ATP model received invalid diameter: "
            f"D = {D:.6e} m"
        )

    if not np.isfinite(Q) or Q <= 0:
        raise ValueError(
            f"ATP model received invalid flow: "
            f"Q = {Q:.6e} m3/s"
        )

    return (
        KD * np.pi * D
        / ((1 - HD) * Q)
    )


def atp_concentration(
    x,
    x0,
    C0,
    S0,
    Q,
    D,
    q
):

    a = alpha(D, S0, q)
    b = beta(D, Q, q)
    g = gamma(D, Q)

    dx = x - x0

    return (
        a
        + b * dx
        + np.exp(-g * dx)
        * (C0 - a)
    )


def atp_profile(
    x,
    x0,
    C0,
    S0,
    Q,
    D,
    q
):

    C = atp_concentration(
        x,
        x0,
        C0,
        S0,
        Q,
        D,
        q
    )

    return x, C
    
    
    print(
    f"{name:16s}",
    f"Q={Q:.6e}",
    f"S0={S0:.6f}",
    f"S1={S[-1]:.6f}",
    f"q={q:.6e}",
    f"C0={C0 * 1000:.6f} uM",
    f"C1={C[-1] * 1000:.6f} uM"
)
