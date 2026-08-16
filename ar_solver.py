from scipy.integrate import solve_ivp

from wall_tension import (
    wall_tension,
    total_tension,
)

from autoregulation import (
    tone_signal,
    target_activation,
    diameter_rate,
    activation_rate,
)


def ar_equations(
    t,
    state,
    P,
    WSS,
    Scr,
    params,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):

    D, A = state

    T = wall_tension(
        P,
        D
    )

    T_total = total_tension(
        D,
        A,
        params["D0"],
        params["C_pass"],
        params["C_pass_prime"],
        params["C_act"],
        params["C_act_prime"],
        params["C_act_double_prime"],
    )

    S_tone = tone_signal(
        T,
        WSS,
        Scr,
        params["C_myo"],
        params["C_tone_myo"],
        params["C_tone_meta"],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
        T_control=params.get("T_control"),
        WSS_control=params.get("WSS_control"),
    )

    A_total = target_activation(
        S_tone
    )

    dDdt = diameter_rate(
        T,
        T_total,
        params["D_control"],
        params["T_control"]
    )

    dAdt = activation_rate(
        A,
        A_total
    )

    return [
        dDdt,
        dAdt
    ]


def solve_ar(
    D_initial,
    A_initial,
    P,
    WSS,
    Scr,
    params,
    t_end=1.0,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):

    return solve_ivp(
        ar_equations,
        (0, t_end),
        [D_initial, A_initial],
        args=(
            P,
            WSS,
            Scr,
            params,
            enable_myo,
            enable_shear,
            enable_meta,
        ),
        rtol=1e-7,
        atol=1e-10,
    )
