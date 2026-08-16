import numpy as np
import pandas as pd

from inputs import INPUT_FILE
from solver import solve_network
from units import um_to_m, cm_to_m, ul_min_to_m3_s, oxygen_demand_to_si
from oxygen import oxygen_consumption, saturation_profile
from ATP import atp_profile, C_INITIAL
from Scr import conducted_response, conducted_profile
from plots import plot_figure2


M0_VALUES = [1.0, 8.28, 20.0]

EFFECTIVE_DIAMETERS_UM = {
    "Artery": 65.2,
    "Vein": 119.1,
}

FULL_PATHWAY = [
    "Artery",
    "Large Arteriole",
    "Small Arteriole",
    "Capillary",
    "Small Venule",
    "Large Venule",
    "Vein",
]

OXYGEN_VESSELS = [
    "Artery",
    "Large Arteriole",
    "Small Arteriole",
    "Capillary",
]

TISSUE_DISTANCE = 19.8e-6
POINTS = 100
S_INITIAL = 0.97


def get_metabolic_diameter(data, name):
    if name in EFFECTIVE_DIAMETERS_UM:
        return um_to_m(EFFECTIVE_DIAMETERS_UM[name])

    row = data[data["compartment"] == name].iloc[0]
    return um_to_m(row["diameter_um"])


def solve_oxygen_atp_segment(
    x, x0, S0, C0, Q, D, L, M0_si, oxygen_active
):
    if not oxygen_active:
        S = np.full(len(x), S0, dtype=float)
        _, C = atp_profile(x, x0, C0, S0, Q, D, 0.0)
        return S, C

    q = oxygen_consumption(M0_si, D, TISSUE_DISTANCE)

    _, S_raw = saturation_profile(
        S0,
        Q,
        D,
        L,
        M0_si,
        TISSUE_DISTANCE,
        len(x),
    )

    if np.all(S_raw >= 0.0):
        S = np.clip(S_raw, 0.0, 1.0)
        _, C = atp_profile(x, x0, C0, S0, Q, D, q)
        return S, C

    first_zero = np.where(S_raw <= 0.0)[0][0]

    if first_zero == 0:
        S = np.zeros(len(x))
        _, C = atp_profile(x, x0, C0, 0.0, Q, D, 0.0)
        return S, C

    i0 = first_zero - 1
    i1 = first_zero

    x_a = x[i0]
    x_b = x[i1]
    S_a = S_raw[i0]
    S_b = S_raw[i1]

    x_zero = x_a + (0.0 - S_a) * (x_b - x_a) / (S_b - S_a)

    x_to_zero = np.append(x[x < x_zero], x_zero)

    _, C_to_zero = atp_profile(
        x_to_zero,
        x0,
        C0,
        S0,
        Q,
        D,
        q,
    )

    C_zero = C_to_zero[-1]
    S = np.clip(S_raw, 0.0, 1.0)
    C = np.empty(len(x))

    before = x <= x_zero
    after = x > x_zero

    if np.any(before):
        _, C_before = atp_profile(
            x[before],
            x0,
            C0,
            S0,
            Q,
            D,
            q,
        )
        C[before] = C_before

    if np.any(after):
        x_after = np.insert(x[after], 0, x_zero)

        _, C_after = atp_profile(
            x_after,
            x_zero,
            C_zero,
            0.0,
            Q,
            D,
            0.0,
        )

        C[after] = C_after[1:]

    return S, C


def solve_metabolic_profile(data, flow_results, M0):
    M0_si = oxygen_demand_to_si(M0)

    x_all = [0.0]
    S_all = [S_INITIAL]
    C_all = [C_INITIAL]

    x0 = 0.0
    S0 = S_INITIAL
    C0 = C_INITIAL

    boundaries = [0.0]
    midpoints = {}
    scr_end = None

    for name in FULL_PATHWAY:
        input_row = data[data["compartment"] == name].iloc[0]
        flow_row = flow_results[
            flow_results["compartment"] == name
        ].iloc[0]

        L = cm_to_m(input_row["length_cm"])
        Q = ul_min_to_m3_s(flow_row["flow_ul_min"])
        D = get_metabolic_diameter(data, name)

        x = np.linspace(x0, x0 + L, POINTS)
        midpoints[name] = x0 + L / 2.0

        S, C = solve_oxygen_atp_segment(
            x=x,
            x0=x0,
            S0=S0,
            C0=C0,
            Q=Q,
            D=D,
            L=L,
            M0_si=M0_si,
            oxygen_active=name in OXYGEN_VESSELS,
        )

        x_all.extend(x[1:])
        S_all.extend(S[1:])
        C_all.extend(C[1:])

        x0 = x[-1]
        S0 = S[-1]
        C0 = C[-1]

        boundaries.append(x0)

        if name == "Large Venule":
            scr_end = x0

    x_all = np.asarray(x_all)
    S_all = np.asarray(S_all)
    C_all = np.asarray(C_all)

    scr_mask = x_all <= scr_end + 1e-12
    x_scr = x_all[scr_mask]
    C_scr = C_all[scr_mask]

    Scr_partial = conducted_profile(x_scr, C_scr)

    Scr_all = np.full(x_all.shape, np.nan)
    Scr_all[scr_mask] = Scr_partial

    Scr_LA = conducted_response(
        x_scr,
        C_scr,
        midpoints["Large Arteriole"],
    )

    Scr_SA = conducted_response(
        x_scr,
        C_scr,
        midpoints["Small Arteriole"],
    )

    return {
        "x": x_all,
        "S": S_all,
        "C": C_all,
        "Scr": Scr_all,
        "Scr_LA": Scr_LA,
        "Scr_SA": Scr_SA,
        "x_LA": midpoints["Large Arteriole"],
        "x_SA": midpoints["Small Arteriole"],
        "x_scr_end": scr_end,
        "boundaries": boundaries,
    }


def print_diagnostics(profile, M0):
    names = [
        "Inlet",
        "Artery",
        "Large Arteriole",
        "Small Arteriole",
        "Capillary",
        "Small Venule",
        "Large Venule",
        "Vein",
    ]

    print(f"\nM0 = {M0}")
    print(
        f"Final S = {profile['S'][-1]:.4f} | "
        f"Final ATP = {profile['C'][-1] * 1000:.4f} uM | "
        f"Scr LA = {profile['Scr_LA'] * 1e5:.4f} | "
        f"Scr SA = {profile['Scr_SA'] * 1e5:.4f}"
    )

    print(
        f"\n{'Location':<20}"
        f"{'x (cm)':>10}"
        f"{'S':>12}"
        f"{'ATP (uM)':>14}"
    )

    for name, xb in zip(names, profile["boundaries"]):
        idx = np.argmin(np.abs(profile["x"] - xb))

        print(
            f"{name:<20}"
            f"{profile['x'][idx] * 100:>10.4f}"
            f"{profile['S'][idx]:>12.6f}"
            f"{profile['C'][idx] * 1000:>14.6f}"
        )


def main():
    data = pd.read_csv(INPUT_FILE)

    flow_results, summary = solve_network(data)

    profiles = {}
    rows = []

    for M0 in M0_VALUES:
        profile = solve_metabolic_profile(
            data,
            flow_results,
            M0,
        )

        profiles[M0] = profile

        for x, S, C, Scr in zip(
            profile["x"],
            profile["S"],
            profile["C"],
            profile["Scr"],
        ):
            rows.append({
                "M0": M0,
                "x_cm": x * 100.0,
                "oxygen_saturation": S,
                "ATP_uM": C * 1000.0,
                "Scr_uM_cm": Scr * 1e5 if np.isfinite(Scr) else np.nan,
            })

    pd.DataFrame(rows).to_csv(
        "results_meta.csv",
        index=False,
    )

    print(f"\nTotal control flow = {summary['Q_total']:.6e} m3/s")
    print("Saved results_meta.csv")

    for M0 in M0_VALUES:
        print_diagnostics(profiles[M0], M0)

    plot_figure2(profiles)


if __name__ == "__main__":
    main()
