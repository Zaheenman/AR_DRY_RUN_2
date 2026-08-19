from scipy.special import expit

from ar_parameters import (
    TAU_D,
    TAU_A,
    C_SHEAR,
    C_META,
)


def _frozen_or_current(
    current_value,
    control_value,
    enabled,
    name,
):

    if enabled:
        return current_value

    if control_value is None:
        raise ValueError(
            f"{name} mechanism is disabled, but no control-state "
            f"value was supplied."
        )

    return control_value


def tone_components(
    T,
    WSS,
    Scr,
    C_myo,
    C_tone_myo,
    C_tone_meta,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
    T_control=None,
    WSS_control=None,
):

    T_used = _frozen_or_current(
        T,
        T_control,
        enable_myo,
        "Myogenic",
    )

    WSS_used = _frozen_or_current(
        WSS,
        WSS_control,
        enable_shear,
        "Shear",
    )

    if enable_meta:
        metabolic = -C_META * Scr
        baseline = C_tone_meta
    else:
        metabolic = 0.0
        baseline = C_tone_myo

    return {
        "myogenic": C_myo * T_used,
        "shear": -C_SHEAR * WSS_used,
        "metabolic": metabolic,
        "baseline": baseline,
    }


def tone_signal(
    T,
    WSS,
    Scr,
    C_myo,
    C_tone_myo,
    C_tone_meta,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
    T_control=None,
    WSS_control=None,
):

    components = tone_components(
        T,
        WSS,
        Scr,
        C_myo,
        C_tone_myo,
        C_tone_meta,
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
        T_control=T_control,
        WSS_control=WSS_control,
    )

    return sum(components.values())


def target_activation(S_tone):

    return expit(S_tone)


def diameter_rate(
    T,
    T_total,
    D_control,
    T_control
):

    return (
        1 / TAU_D
        * D_control / T_control
        * (T - T_total)
    )


def activation_rate(
    A,
    A_total
):

    return (
        1 / TAU_A
        * (A_total - A)
    )
