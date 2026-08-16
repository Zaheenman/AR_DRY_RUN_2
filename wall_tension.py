import numpy as np


def wall_tension(P, D):

    return P * D / 2


def passive_tension(
    D,
    D0,
    C_pass,
    C_pass_prime
):

    return C_pass * np.exp(
        C_pass_prime
        * (D / D0 - 1)
    )


def active_tension_max(
    D,
    D0,
    C_act,
    C_act_prime,
    C_act_double_prime
):

    return C_act * np.exp(
        -(
            (D / D0 - C_act_prime)
            / C_act_double_prime
        )**2
    )


def total_tension(
    D,
    A,
    D0,
    C_pass,
    C_pass_prime,
    C_act,
    C_act_prime,
    C_act_double_prime
):

    T_pass = passive_tension(
        D,
        D0,
        C_pass,
        C_pass_prime
    )

    T_active = active_tension_max(
        D,
        D0,
        C_act,
        C_act_prime,
        C_act_double_prime
    )

    return (
        T_pass
        + A * T_active
    )
