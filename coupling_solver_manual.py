import numpy as np
import pandas as pd

from oxygen import oxygen_consumption, saturation_profile
from ATP import atp_profile, C_INITIAL
from Scr import conducted_response
from ar_solver import solve_ar
from ar_parameters import AR_PARAMS
from autoregulation import tone_components, tone_signal, target_activation
from wall_tension import wall_tension, total_tension
from flow import velocity, wall_shear_stress
from units import (
    um_to_m,
    cm_to_m,
    cp_to_pa_s,
    ul_min_to_m3_s,
    mmhg_to_pa,
    oxygen_demand_to_si,
    m_s_to_cm_s,
)


ACTIVE_VESSELS = [
    "Large Arteriole",
    "Small Arteriole",
]


METABOLIC_VESSELS = [
    "Large Arteriole",
    "Small Arteriole",
    "Capillary",
    "Small Venule",
    "Large Venule",
]


OXYGEN_VESSELS = [
    "Large Arteriole",
    "Small Arteriole",
    "Capillary",
]


TISSUE_DISTANCE = 18.8e-6
S_INITIAL = 0.97


def _input_row(data, vessel):
    row = data.loc[
        data["compartment"] == vessel
    ]

    if row.empty:
        raise ValueError(
            f"Compartment '{vessel}' not found in input data."
        )

    return row.iloc[0]


def _manual_value(manual_inputs, vessel, key):
    if vessel not in manual_inputs:
        raise ValueError(
            f"Manual haemodynamic input missing for '{vessel}'."
        )

    if key not in manual_inputs[vessel]:
        raise ValueError(
            f"Manual input '{key}' missing for '{vessel}'."
        )

    value = manual_inputs[vessel][key]

    if value is None or not np.isfinite(value):
        raise ValueError(
            f"Invalid manual {key} for '{vessel}': {value}"
        )

    return float(value)


def build_manual_results(
    data,
    manual_inputs,
):
    rows = []

    for vessel, values in manual_inputs.items():
        input_row = _input_row(data, vessel)

        pressure_mmhg = _manual_value(
            manual_inputs,
            vessel,
            "pressure_mmhg",
        )

        flow_ul_min = _manual_value(
            manual_inputs,
            vessel,
            "flow_ul_min",
        )

        diameter_um = float(
            input_row["diameter_um"]
        )

        viscosity_cp = float(
            input_row["viscosity_cp"]
        )

        D = um_to_m(diameter_um)
        mu = cp_to_pa_s(viscosity_cp)
        Q = ul_min_to_m3_s(flow_ul_min)

        rows.append(
            {
                "compartment": vessel,
                "diameter_um": diameter_um,
                "pressure_mmhg": pressure_mmhg,
                "flow_ul_min": flow_ul_min,
                "velocity_cm_s": m_s_to_cm_s(
                    velocity(Q, D)
                ),
                "wss_pa": wall_shear_stress(
                    mu,
                    Q,
                    D,
                ),
            }
        )

    return pd.DataFrame(rows)


def metabolic_profiles(
    data,
    manual_results,
    M0,
    points=100,
):
    M0_si = oxygen_demand_to_si(M0)

    artery_length = cm_to_m(
        _input_row(
            data,
            "Artery",
        )["length_cm"]
    )

    x_all = [
        0.0,
        artery_length,
    ]

    S_all = [
        S_INITIAL,
        S_INITIAL,
    ]

    C_all = [
        C_INITIAL,
        C_INITIAL,
    ]

    x0 = artery_length
    S0 = S_INITIAL
    C0 = C_INITIAL

    midpoints = {}

    for name in METABOLIC_VESSELS:
        input_row = _input_row(
            data,
            name,
        )

        flow_row = manual_results.loc[
            manual_results["compartment"] == name
        ]

        if flow_row.empty:
            raise ValueError(
                f"Manual Q/P data missing for metabolic vessel '{name}'."
            )

        flow_row = flow_row.iloc[0]

        D = um_to_m(
            input_row["diameter_um"]
        )

        L = cm_to_m(
            input_row["length_cm"]
        )

        Q = ul_min_to_m3_s(
            flow_row["flow_ul_min"]
        )

        if Q <= 0.0:
            raise ValueError(
                f"{name}: manual flow must be positive."
            )

        x = np.linspace(
            x0,
            x0 + L,
            points,
        )

        midpoints[name] = (
            x0 + L / 2
        )

        if name in OXYGEN_VESSELS:
            q = oxygen_consumption(
                M0_si,
                D,
                TISSUE_DISTANCE,
            )

            _, S_raw = saturation_profile(
                S0,
                Q,
                D,
                L,
                M0_si,
                TISSUE_DISTANCE,
                points,
            )

            S = np.clip(
                S_raw,
                0.0,
                1.0,
            )

        else:
            q = 0.0

            S_raw = np.full(
                points,
                S0,
            )

            S = S_raw.copy()

        _, C = atp_profile(
            x,
            x0,
            C0,
            S0,
            Q,
            D,
            q,
        )

        x_all.extend(
            x[1:]
        )

        S_all.extend(
            S[1:]
        )

        C_all.extend(
            C[1:]
        )

        x0 = x[-1]
        S0 = S[-1]
        C0 = C[-1]

    x_all = np.array(
        x_all
    )

    S_all = np.array(
        S_all
    )

    C_all = np.array(
        C_all
    )

    Scr_LA = conducted_response(
        x_all,
        C_all,
        midpoints["Large Arteriole"],
    )

    Scr_SA = conducted_response(
        x_all,
        C_all,
        midpoints["Small Arteriole"],
    )

    return {
        "x": x_all,
        "S": S_all,
        "C": C_all,
        "Scr_LA": Scr_LA,
        "Scr_SA": Scr_SA,
        "S_out": S_all[-1],
        "ATP_out": C_all[-1],
        "ATP_out_uM": C_all[-1] * 1000.0,
    }


def control_tension(
    manual_results,
    vessel,
):
    row = manual_results.loc[
        manual_results["compartment"] == vessel
    ].iloc[0]

    P = mmhg_to_pa(
        row["pressure_mmhg"]
    )

    D = um_to_m(
        row["diameter_um"]
    )

    return wall_tension(
        P,
        D,
    )


def build_parameters(
    control_results,
    control_metabolic,
    A_control=0.6,
    myo_factor=1.0,
):
    parameters = {
        name: AR_PARAMS[name].copy()
        for name in ACTIVE_VESSELS
    }

    S_control = np.log(
        A_control / (1.0 - A_control)
    )

    for name in ACTIVE_VESSELS:
        row = control_results.loc[
            control_results["compartment"] == name
        ].iloc[0]

        T_control = control_tension(
            control_results,
            name,
        )

        WSS_control = row["wss_pa"]

        Scr_control = (
            control_metabolic["Scr_LA"]
            if name == "Large Arteriole"
            else control_metabolic["Scr_SA"]
        )

        parameters[name]["C_myo"] *= myo_factor

        parameters[name]["D_control"] = (
            row["diameter_um"] * 1e-6
        )

        parameters[name]["A_control"] = A_control
        parameters[name]["T_control"] = T_control
        parameters[name]["WSS_control"] = WSS_control

        parameters[name]["C_tone_myo"] = (
            S_control
            - parameters[name]["C_myo"] * T_control
            + C_SHEAR * WSS_control
        )

        parameters[name]["C_tone_meta"] = (
            S_control
            - parameters[name]["C_myo"] * T_control
            + C_SHEAR * WSS_control
            + C_META * Scr_control
        )

    return parameters


def _ensure_control_stimuli(
    control_results,
    parameters,
):
    for name in ACTIVE_VESSELS:
        if parameters[name].get(
            "T_control"
        ) is None:
            parameters[name]["T_control"] = (
                control_tension(
                    control_results,
                    name,
                )
            )

        if parameters[name].get(
            "WSS_control"
        ) is None:
            row = control_results.loc[
                control_results["compartment"] == name
            ].iloc[0]

            parameters[name]["WSS_control"] = (
                row["wss_pa"]
            )


def advance_vessel(
    manual_results,
    vessel,
    activation,
    Scr,
    params,
    dt,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):
    row = manual_results.loc[
        manual_results["compartment"] == vessel
    ].iloc[0]

    D = um_to_m(
        row["diameter_um"]
    )

    P = mmhg_to_pa(
        row["pressure_mmhg"]
    )

    WSS = row["wss_pa"]

    solution = solve_ar(
        D,
        activation,
        P,
        WSS,
        Scr,
        params,
        t_end=dt,
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
    )

    return (
        solution.y[0, -1],
        solution.y[1, -1],
    )


def final_vessel_state(
    manual_results,
    vessel,
    activation,
    Scr,
    params,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):
    row = manual_results.loc[
        manual_results["compartment"] == vessel
    ].iloc[0]

    D = um_to_m(
        row["diameter_um"]
    )

    P = mmhg_to_pa(
        row["pressure_mmhg"]
    )

    WSS = row["wss_pa"]

    T = wall_tension(
        P,
        D,
    )

    T_total = total_tension(
        D,
        activation,
        params["D0"],
        params["C_pass"],
        params["C_pass_prime"],
        params["C_act"],
        params["C_act_prime"],
        params["C_act_double_prime"],
    )

    components = tone_components(
        T,
        WSS,
        Scr,
        params["C_myo"],
        params["C_tone_myo"],
        params["C_tone_meta"],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
        T_control=params.get(
            "T_control"
        ),
        WSS_control=params.get(
            "WSS_control"
        ),
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
        T_control=params.get(
            "T_control"
        ),
        WSS_control=params.get(
            "WSS_control"
        ),
    )

    return {
        "pressure_mmhg": row["pressure_mmhg"],
        "flow_ul_min": row["flow_ul_min"],
        "wss_pa": WSS,
        "T": T,
        "T_total": T_total,
        "S_tone": S_tone,
        "A_target": target_activation(
            S_tone
        ),
        "myogenic": components["myogenic"],
        "shear": components["shear"],
        "metabolic": components["metabolic"],
        "baseline": components["baseline"],
    }


def _manual_summary(
    manual_inputs,
    Q_total_ul_min,
):
    return {
        "Q_total": ul_min_to_m3_s(
            Q_total_ul_min
        ),
        "Q_total_ul_min": float(
            Q_total_ul_min
        ),
        "manual_inputs": manual_inputs,
    }


def solve_coupled_ar(
    data,
    M0,
    manual_inputs,
    Q_total_ul_min,
    dt=1.0,
    max_iterations=1000,
    diameter_tolerance=1e-9,
    activation_tolerance=1e-6,
    control_results=None,
    parameters=None,
    initial_activations=None,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):
    data = data.copy()

    data["compartment"] = (
        data["compartment"]
        .astype(str)
        .str.strip()
    )

    if initial_activations is None:
        activations = {
            name: AR_PARAMS[name][
                "A_control"
            ]
            for name in ACTIVE_VESSELS
        }
    else:
        activations = (
            initial_activations.copy()
        )

    if control_results is None:
        control_results = (
            build_manual_results(
                data,
                manual_inputs,
            )
        )

    if parameters is None:
        parameters = build_parameters(
            control_results
        )

    _ensure_control_stimuli(
        control_results,
        parameters,
    )

    converged = False
    diameter_error = np.inf
    activation_error = np.inf

    for iteration in range(
        max_iterations
    ):
        manual_results = (
            build_manual_results(
                data,
                manual_inputs,
            )
        )

        metabolic = (
            metabolic_profiles(
                data,
                manual_results,
                M0,
            )
        )

        old_diameters = {}
        new_diameters = {}
        new_activations = {}

        for name in ACTIVE_VESSELS:
            old_diameters[name] = (
                um_to_m(
                    data.loc[
                        data[
                            "compartment"
                        ] == name,
                        "diameter_um",
                    ].iloc[0]
                )
            )

            Scr = (
                metabolic["Scr_LA"]
                if name == "Large Arteriole"
                else metabolic["Scr_SA"]
            )

            D_new, A_new = (
                advance_vessel(
                    manual_results,
                    name,
                    activations[name],
                    Scr,
                    parameters[name],
                    dt,
                    enable_myo=enable_myo,
                    enable_shear=enable_shear,
                    enable_meta=enable_meta,
                )
            )

            if (
                not np.isfinite(D_new)
                or D_new <= 0.0
            ):
                raise ValueError(
                    f"{name}: autoregulation produced "
                    f"invalid D_new={D_new:.6e} m"
                )

            new_diameters[name] = (
                D_new
            )

            new_activations[name] = (
                A_new
            )

        diameter_error = max(
            abs(
                new_diameters[name]
                - old_diameters[name]
            )
            for name in ACTIVE_VESSELS
        )

        activation_error = max(
            abs(
                new_activations[name]
                - activations[name]
            )
            for name in ACTIVE_VESSELS
        )

        for name in ACTIVE_VESSELS:
            data.loc[
                data[
                    "compartment"
                ] == name,
                "diameter_um",
            ] = (
                new_diameters[name]
                * 1e6
            )

            activations[name] = (
                new_activations[name]
            )

        if (
            diameter_error
            < diameter_tolerance
            and activation_error
            < activation_tolerance
        ):
            converged = True
            break

    manual_results = (
        build_manual_results(
            data,
            manual_inputs,
        )
    )

    metabolic = metabolic_profiles(
        data,
        manual_results,
        M0,
    )

    LA_state = final_vessel_state(
        manual_results,
        "Large Arteriole",
        activations[
            "Large Arteriole"
        ],
        metabolic["Scr_LA"],
        parameters[
            "Large Arteriole"
        ],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
    )

    SA_state = final_vessel_state(
        manual_results,
        "Small Arteriole",
        activations[
            "Small Arteriole"
        ],
        metabolic["Scr_SA"],
        parameters[
            "Small Arteriole"
        ],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
    )

    summary = _manual_summary(
        manual_inputs,
        Q_total_ul_min,
    )

    return {
        "data": data,
        "flow_results": manual_results,
        "manual_results": manual_results,
        "summary": summary,
        "metabolic": metabolic,
        "diameter_LA_um": data.loc[
            data[
                "compartment"
            ] == "Large Arteriole",
            "diameter_um",
        ].iloc[0],
        "diameter_SA_um": data.loc[
            data[
                "compartment"
            ] == "Small Arteriole",
            "diameter_um",
        ].iloc[0],
        "activation_LA": activations[
            "Large Arteriole"
        ],
        "activation_SA": activations[
            "Small Arteriole"
        ],
        "Scr_LA": metabolic[
            "Scr_LA"
        ],
        "Scr_SA": metabolic[
            "Scr_SA"
        ],
        "ATP_out": metabolic[
            "ATP_out"
        ],
        "ATP_out_uM": metabolic[
            "ATP_out_uM"
        ],
        "S_out": metabolic[
            "S_out"
        ],
        "LA_state": LA_state,
        "SA_state": SA_state,
        "iterations": iteration + 1,
        "converged": converged,
        "diameter_error": diameter_error,
        "activation_error": activation_error,
    }


def solve_figure4(
    data,
    M0_values,
    manual_inputs,
    Q_total_ul_min,
    **solver_options,
):
    control_results = (
        build_manual_results(
            data,
            manual_inputs,
        )
    )

    parameters = build_parameters(
        control_results
    )

    rows = []

    for M0 in M0_values:
        result = solve_coupled_ar(
            data,
            M0,
            manual_inputs,
            Q_total_ul_min,
            control_results=control_results,
            parameters=parameters,
            **solver_options,
        )

        LA = result[
            "LA_state"
        ]

        SA = result[
            "SA_state"
        ]

        rows.append(
            {
                "M0": M0,
                "diameter_LA_um": result[
                    "diameter_LA_um"
                ],
                "diameter_SA_um": result[
                    "diameter_SA_um"
                ],
                "activation_LA": result[
                    "activation_LA"
                ],
                "activation_SA": result[
                    "activation_SA"
                ],
                "Scr_LA_uM_cm": (
                    result["Scr_LA"]
                    * 1e5
                ),
                "Scr_SA_uM_cm": (
                    result["Scr_SA"]
                    * 1e5
                ),
                "flow_ul_min": (
                    Q_total_ul_min
                ),
                "ATP_out_uM": result[
                    "ATP_out_uM"
                ],
                "S_out": result[
                    "S_out"
                ],
                "LA_myogenic": LA[
                    "myogenic"
                ],
                "LA_shear": LA[
                    "shear"
                ],
                "LA_metabolic": LA[
                    "metabolic"
                ],
                "LA_baseline": LA[
                    "baseline"
                ],
                "LA_S_tone": LA[
                    "S_tone"
                ],
                "SA_myogenic": SA[
                    "myogenic"
                ],
                "SA_shear": SA[
                    "shear"
                ],
                "SA_metabolic": SA[
                    "metabolic"
                ],
                "SA_baseline": SA[
                    "baseline"
                ],
                "SA_S_tone": SA[
                    "S_tone"
                ],
                "iterations": result[
                    "iterations"
                ],
                "converged": result[
                    "converged"
                ],
            }
        )

    return pd.DataFrame(
        rows
    )
