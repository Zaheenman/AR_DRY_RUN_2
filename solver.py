import pandas as pd

from units import *
from flow import *


def solve_network(data):

    data = data.copy()

    p_in = mmhg_to_pa(data["p_in_mmhg"].dropna().iloc[0])
    p_out = mmhg_to_pa(data["p_out_mmhg"].dropna().iloc[0])

    data["D_m"] = data["diameter_um"].apply(
        lambda x: um_to_m(x) if pd.notna(x) else None
    )

    data["L_m"] = data["length_cm"].apply(cm_to_m)
    data["mu_Pa_s"] = data["viscosity_cp"].apply(cp_to_pa_s)

    resistances = []

    for _, row in data.iterrows():

        if pd.notna(row["diameter_um"]) and pd.notna(row["n_vessels"]):

            R = resistance(
                row["L_m"],
                row["mu_Pa_s"],
                row["D_m"],
                row["n_vessels"]
            )

        else:

            R = resistance_table_to_si(
                row["resistance"]
            )

        resistances.append(R)

    data["resistance_si"] = resistances

    R_total = total_resistance(
        data["resistance_si"]
    )

    Q_total = total_flow(
        p_in,
        p_out,
        R_total
    )

    current_pressure = p_in

    pressures = []
    flows = []
    velocities = []
    wss = []

    for _, row in data.iterrows():

        dP = pressure_drop(
            Q_total,
            row["resistance_si"]
        )

        p_end = current_pressure - dP
        p_mean = mean_pressure(
            current_pressure,
            p_end
        )

        pressures.append(
            pa_to_mmhg(p_mean)
        )

        if pd.notna(row["diameter_um"]) and pd.notna(row["n_vessels"]):

            Q = vessel_flow(
                Q_total,
                row["n_vessels"]
            )

            v = velocity(
                Q,
                row["D_m"]
            )

            tau = wall_shear_stress(
                row["mu_Pa_s"],
                Q,
                row["D_m"]
            )

            flows.append(
                m3_s_to_ul_min(Q)
            )

            velocities.append(
                m_s_to_cm_s(v)
            )

            wss.append(tau)

        else:

            flows.append(
                m3_s_to_ul_min(Q_total)
            )

            velocities.append(None)
            wss.append(None)

        current_pressure = p_end

    data["pressure_mmhg"] = pressures
    data["flow_ul_min"] = flows
    data["velocity_cm_s"] = velocities
    data["wss_pa"] = wss

    results = data[
        [
            "compartment",
            "diameter_um",
            "resistance_si",
            "pressure_mmhg",
            "flow_ul_min",
            "velocity_cm_s",
            "wss_pa",
        ]
    ]

    summary = {
        "R_total": R_total,
        "Q_total": Q_total,
        "p_in": p_in,
        "p_out": current_pressure,
    }

    return results, summary
