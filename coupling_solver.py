import numpy as np
import pandas as pd

from solver import solve_network
from oxygen import oxygen_consumption, saturation_profile
from ATP import atp_profile, C_INITIAL
from Scr import conducted_response
from ar_solver import solve_ar
from ar_parameters import AR_PARAMS
from autoregulation import tone_components, tone_signal, target_activation
from wall_tension import wall_tension, total_tension
from units import (
    um_to_m,
    cm_to_m,
    ul_min_to_m3_s,
    mmhg_to_pa,
    oxygen_demand_to_si,
)


ACTIVE_VESSELS = [
    "Large Arteriole",
    "Small Arteriole"
]

METABOLIC_VESSELS = [
    "Large Arteriole",
    "Small Arteriole",
    "Capillary",
    "Small Venule",
    "Large Venule"
]

OXYGEN_VESSELS = [
    "Large Arteriole",
    "Small Arteriole",
    "Capillary"
]

TISSUE_DISTANCE = 18.8e-6
S_INITIAL = 0.97


def metabolic_profiles(
    data,
    flow_results,
    M0,
    points=100
):

    M0_si = oxygen_demand_to_si(M0)

    artery_length = cm_to_m(
        data.loc[
            data["compartment"] == "Artery",
            "length_cm"
        ].iloc[0]
    )

    x_all = [0.0, artery_length]
    S_all = [S_INITIAL, S_INITIAL]
    C_all = [C_INITIAL, C_INITIAL]

    x0 = artery_length
    S0 = S_INITIAL
    C0 = C_INITIAL

    midpoints = {}

    for name in METABOLIC_VESSELS:

        input_row = data[
            data["compartment"] == name
        ].iloc[0]

        flow_row = flow_results[
            flow_results["compartment"] == name
        ].iloc[0]

        D = um_to_m(
            input_row["diameter_um"]
        )

        L = cm_to_m(
            input_row["length_cm"]
        )

        Q = ul_min_to_m3_s(
            flow_row["flow_ul_min"]
        )

        x = np.linspace(
            x0,
            x0 + L,
            points
        )

        midpoints[name] = x0 + L / 2

        if name in OXYGEN_VESSELS:

            q = oxygen_consumption(
                M0_si,
                D,
                TISSUE_DISTANCE
            )

            _, S = saturation_profile(
                S0,
                Q,
                D,
                L,
                M0_si,
                TISSUE_DISTANCE,
                points
            )

            S = np.clip(
                S,
                0.0,
                1.0
            )

        else:

            q = 0.0

            S = np.full(
                points,
                S0
            )

        _, C = atp_profile(
            x,
            x0,
            C0,
            S0,
            Q,
            D,
            q
        )

        x_all.extend(x[1:])
        S_all.extend(S[1:])
        C_all.extend(C[1:])

        x0 = x[-1]
        S0 = S[-1]
        C0 = C[-1]

    x_all = np.array(x_all)
    S_all = np.array(S_all)
    C_all = np.array(C_all)

    Scr_LA = conducted_response(
        x_all,
        C_all,
        midpoints["Large Arteriole"]
    )

    Scr_SA = conducted_response(
        x_all,
        C_all,
        midpoints["Small Arteriole"]
    )

    return {
        "x": x_all,
        "S": S_all,
        "C": C_all,
        "Scr_LA": Scr_LA,
        "Scr_SA": Scr_SA,
    }


def control_tension(
    flow_results,
    vessel
):

    row = flow_results[
        flow_results["compartment"] == vessel
    ].iloc[0]

    P = mmhg_to_pa(
        row["pressure_mmhg"]
    )

    D = um_to_m(
        row["diameter_um"]
    )

    return wall_tension(P, D)


def build_parameters(
    control_results
):
    """
    Build the AR parameter dictionaries from the haemodynamic control state.

    The paper uses two VSM constants:
      C_tone_myo  -> Eq. 6, without the metabolic term
      C_tone_meta -> Eq. 7, with the metabolic term

    T_control and WSS_control are stored so that myogenic and shear
    responses can be deactivated by holding their stimuli at control.
    """

    parameters = {
        name: AR_PARAMS[name].copy()
        for name in ACTIVE_VESSELS
    }

    for name in ACTIVE_VESSELS:

        row = control_results[
            control_results["compartment"] == name
        ].iloc[0]

        parameters[name]["T_control"] = control_tension(
            control_results,
            name
        )

        parameters[name]["WSS_control"] = row["wss_pa"]

    return parameters


def _ensure_control_stimuli(
    data,
    control_results,
    parameters,
    M0
):
    """
    Ensure the frozen control stimuli required for deactivated
    myogenic/shear responses are available.

    data and M0 remain in the signature to preserve the existing
    coupling-solver call structure.
    """

    for name in ACTIVE_VESSELS:

        if parameters[name].get("T_control") is None:
            parameters[name]["T_control"] = control_tension(
                control_results,
                name
            )

        if parameters[name].get("WSS_control") is None:

            row = control_results[
                control_results["compartment"] == name
            ].iloc[0]

            parameters[name]["WSS_control"] = row["wss_pa"]


def advance_vessel(
    flow_results,
    vessel,
    activation,
    Scr,
    params,
    dt,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):

    row = flow_results[
        flow_results["compartment"] == vessel
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
        solution.y[1, -1]
    )


def final_vessel_state(
    flow_results,
    vessel,
    activation,
    Scr,
    params,
    enable_myo=True,
    enable_shear=True,
    enable_meta=True,
):

    row = flow_results[
        flow_results["compartment"] == vessel
    ].iloc[0]

    D = um_to_m(row["diameter_um"])
    P = mmhg_to_pa(row["pressure_mmhg"])
    WSS = row["wss_pa"]

    T = wall_tension(P, D)

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
        T_control=params.get("T_control"),
        WSS_control=params.get("WSS_control"),
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

    return {
        "pressure_mmhg": row["pressure_mmhg"],
        "wss_pa": WSS,
        "T": T,
        "T_total": T_total,
        "S_tone": S_tone,
        "A_target": target_activation(S_tone),
        "myogenic": components["myogenic"],
        "shear": components["shear"],
        "metabolic": components["metabolic"],
        "baseline": components["baseline"],
    }


def solve_coupled_ar(
    data,
    M0,
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

    if initial_activations is None:
        activations = {
            name: AR_PARAMS[name]["A_control"]
            for name in ACTIVE_VESSELS
        }
    else:
        activations = initial_activations.copy()

    if control_results is None:
        control_results, _ = solve_network(data)

    if parameters is None:
        parameters = build_parameters(control_results)

    # Important:
    # Disabled mechanisms are frozen at their control stimulus, not removed.
    # Populate those control stimuli before the first AR update.
    _ensure_control_stimuli(
        data,
        control_results,
        parameters,
        M0
    )

    converged = False
    diameter_error = np.inf
    activation_error = np.inf

    for iteration in range(max_iterations):

        flow_results, summary = solve_network(data)

        # Keep calculating the metabolic state even when metabolic regulation
        # is disabled so ATP/Scr remain available for diagnostics. When
        # enable_meta=False, the tone equation omits the metabolic term and uses
        # the non-metabolic C_tone value from Eq. 6.
        metabolic = metabolic_profiles(
            data,
            flow_results,
            M0
        )

        old_diameters = {}
        new_diameters = {}
        new_activations = {}

        for name in ACTIVE_VESSELS:

            old_diameters[name] = um_to_m(
                data.loc[
                    data["compartment"] == name,
                    "diameter_um"
                ].iloc[0]
            )

            Scr = (
                metabolic["Scr_LA"]
                if name == "Large Arteriole"
                else metabolic["Scr_SA"]
            )

            D_new, A_new = advance_vessel(
                flow_results,
                name,
                activations[name],
                Scr,
                parameters[name],
                dt,
                enable_myo=enable_myo,
                enable_shear=enable_shear,
                enable_meta=enable_meta,
            )

            if not np.isfinite(D_new) or D_new <= 0.0:
                raise ValueError(
                    f"{name}: autoregulation produced an invalid diameter "
                    f"before the next flow/metabolic solve: "
                    f"D_new={D_new:.6e} m, "
                    f"D_old={old_diameters[name]:.6e} m, "
                    f"A_old={activations[name]:.6e}, "
                    f"A_new={A_new:.6e}, "
                    f"enable_myo={enable_myo}, "
                    f"enable_shear={enable_shear}, "
                    f"enable_meta={enable_meta}"
                )

            new_diameters[name] = D_new
            new_activations[name] = A_new

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
                data["compartment"] == name,
                "diameter_um"
            ] = new_diameters[name] * 1e6

            activations[name] = new_activations[name]

        if (
            diameter_error < diameter_tolerance
            and activation_error < activation_tolerance
        ):
            converged = True
            break

    flow_results, summary = solve_network(data)

    metabolic = metabolic_profiles(
        data,
        flow_results,
        M0
    )

    LA_state = final_vessel_state(
        flow_results,
        "Large Arteriole",
        activations["Large Arteriole"],
        metabolic["Scr_LA"],
        parameters["Large Arteriole"],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
    )

    SA_state = final_vessel_state(
        flow_results,
        "Small Arteriole",
        activations["Small Arteriole"],
        metabolic["Scr_SA"],
        parameters["Small Arteriole"],
        enable_myo=enable_myo,
        enable_shear=enable_shear,
        enable_meta=enable_meta,
    )

    return {
        "data": data,
        "flow_results": flow_results,
        "summary": summary,
        "diameter_LA_um": data.loc[
            data["compartment"] == "Large Arteriole",
            "diameter_um"
        ].iloc[0],
        "diameter_SA_um": data.loc[
            data["compartment"] == "Small Arteriole",
            "diameter_um"
        ].iloc[0],
        "activation_LA": activations["Large Arteriole"],
        "activation_SA": activations["Small Arteriole"],
        "Scr_LA": metabolic["Scr_LA"],
        "Scr_SA": metabolic["Scr_SA"],
        "LA_state": LA_state,
        "SA_state": SA_state,
        "iterations": iteration + 1,
        "converged": converged,
    }


def solve_figure4(
    data,
    M0_values,
    **solver_options
):

    control_results, _ = solve_network(data)
    parameters = build_parameters(control_results)

    rows = []

    for M0 in M0_values:

        result = solve_coupled_ar(
            data,
            M0,
            control_results=control_results,
            parameters=parameters,
            **solver_options
        )

        LA = result["LA_state"]
        SA = result["SA_state"]

        rows.append({
            "M0": M0,
            "diameter_LA_um": result["diameter_LA_um"],
            "diameter_SA_um": result["diameter_SA_um"],
            "activation_LA": result["activation_LA"],
            "activation_SA": result["activation_SA"],
            "Scr_LA_uM_cm": result["Scr_LA"] * 1e5,
            "Scr_SA_uM_cm": result["Scr_SA"] * 1e5,
            "flow_ul_min": result["summary"]["Q_total"] * 60 * 1e9,
            "LA_myogenic": LA["myogenic"],
            "LA_shear": LA["shear"],
            "LA_metabolic": LA["metabolic"],
            "LA_baseline": LA["baseline"],
            "LA_S_tone": LA["S_tone"],
            "SA_myogenic": SA["myogenic"],
            "SA_shear": SA["shear"],
            "SA_metabolic": SA["metabolic"],
            "SA_baseline": SA["baseline"],
            "SA_S_tone": SA["S_tone"],
            "iterations": result["iterations"],
            "converged": result["converged"],
        })

    return pd.DataFrame(rows)
