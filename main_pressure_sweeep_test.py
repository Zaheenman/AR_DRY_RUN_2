import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from inputs import INPUT_FILE
from solver import solve_network
from coupling_solver import solve_coupled_ar, build_parameters


# ============================================================
# SETTINGS
# ============================================================

M0 = 8.28

CONTROL_PRESSURE = 100.0

MIN_PRESSURE = 20.0
MAX_PRESSURE = 200.0
PRESSURE_STEP = 10.0

DT = 0.1
MAX_ITERATIONS = 5000


# ============================================================
# MECHANISM STATES
# ============================================================

MECHANISM_STATES = {

    "All": {
        "enable_myo": True,
        "enable_shear": True,
        "enable_meta": True,
    },

    "Myo + Shear": {
        "enable_myo": True,
        "enable_shear": True,
        "enable_meta": False,
    },

    "Myo": {
        "enable_myo": True,
        "enable_shear": False,
        "enable_meta": False,
    },
}


# ============================================================
# SET PRESSURE
# ============================================================

def set_inlet_pressure(data, pressure):

    data = data.copy()

    data.loc[
        data["compartment"] == "Artery",
        "p_in_mmhg"
    ] = pressure

    return data


# ============================================================
# EXTRACT RESULTS
# ============================================================

def extract_row(
    pressure,
    result,
    mechanism,
):

    Q_ul_min = (
        result["summary"]["Q_total"]
        * 60.0
        * 1e9
    )

    LA = result["LA_state"]
    SA = result["SA_state"]

    return {

        "mechanism": mechanism,

        "pressure_mmhg": pressure,

        "flow_ul_min": Q_ul_min,

        "diameter_LA_um":
            result["diameter_LA_um"],

        "diameter_SA_um":
            result["diameter_SA_um"],

        "activation_LA":
            result["activation_LA"],

        "activation_SA":
            result["activation_SA"],

        "pressure_LA_mmhg":
            LA["pressure_mmhg"],

        "pressure_SA_mmhg":
            SA["pressure_mmhg"],

        "wss_LA_pa":
            LA["wss_pa"],

        "wss_SA_pa":
            SA["wss_pa"],

        "Scr_LA":
            result["Scr_LA"],

        "Scr_SA":
            result["Scr_SA"],

        "LA_myogenic":
            LA["myogenic"],

        "LA_shear":
            LA["shear"],

        "LA_metabolic":
            LA["metabolic"],

        "LA_S_tone":
            LA["S_tone"],

        "SA_myogenic":
            SA["myogenic"],

        "SA_shear":
            SA["shear"],

        "SA_metabolic":
            SA["metabolic"],

        "SA_S_tone":
            SA["S_tone"],

        "iterations":
            result["iterations"],

        "converged":
            result["converged"],
    }


# ============================================================
# SOLVE ONE PRESSURE
# ============================================================

def solve_pressure(
    pressure,
    data_start,
    activations_start,
    control_results,
    parameters,
    settings,
):

    data_test = set_inlet_pressure(
        data_start,
        pressure
    )

    return solve_coupled_ar(
        data_test,
        M0,

        control_results=control_results,
        parameters=parameters,

        initial_activations=activations_start,

        enable_myo=settings["enable_myo"],
        enable_shear=settings["enable_shear"],
        enable_meta=settings["enable_meta"],

        dt=DT,
        max_iterations=MAX_ITERATIONS,
    )


# ============================================================
# PRINT RESULT
# ============================================================

def print_result(
    pressure,
    result,
):

    Q_ul_min = (
        result["summary"]["Q_total"]
        * 60.0
        * 1e9
    )

    print(
        f"P={pressure:.0f} mmHg | "
        f"Q={Q_ul_min:.6f} uL/min | "
        f"LA D={result['diameter_LA_um']:.3f} um | "
        f"SA D={result['diameter_SA_um']:.3f} um | "
        f"LA A={result['activation_LA']:.4f} | "
        f"SA A={result['activation_SA']:.4f} | "
        f"converged={result['converged']}"
    )


# ============================================================
# RUN ONE COMPLETELY INDEPENDENT MECHANISM STATE
# ============================================================

def run_pressure_sweep(
    mechanism_name,
    settings,
):

    print()
    print("================================================")
    print(f" MECHANISM: {mechanism_name}")
    print("================================================")

    print(
        f"Myogenic  = {settings['enable_myo']}\n"
        f"Shear     = {settings['enable_shear']}\n"
        f"Metabolic = {settings['enable_meta']}"
    )

    # ========================================================
    # HARD RESET
    #
    # Reload the ORIGINAL input file for every mechanism.
    #
    # Do not use data/results from the previous mechanism.
    # ========================================================

    data_control = pd.read_csv(
        INPUT_FILE
    )

    data_control = set_inlet_pressure(
        data_control,
        CONTROL_PRESSURE
    )

    # ========================================================
    # REBUILD CONTROL HEMODYNAMICS
    # ========================================================

    control_results, control_summary = (
        solve_network(data_control)
    )

    # ========================================================
    # REBUILD PARAMETERS FROM ORIGINAL CONTROL STATE
    # ========================================================

    parameters = build_parameters(
        control_results
    )

    print()
    print("Fresh control state")

    print(
        f"Control inlet pressure = "
        f"{CONTROL_PRESSURE:.1f} mmHg"
    )

    print(
        f"Control flow = "
        f"{control_summary['Q_total'] * 60e9:.6f} "
        f"uL/min"
    )

    # ========================================================
    # START THIS MECHANISM FROM 100 mmHg
    # ========================================================

    control_ar = solve_coupled_ar(
        data_control,
        M0,

        control_results=control_results,
        parameters=parameters,

        # IMPORTANT:
        # Do not pass activation from another run.
        # This starts from the default control activation.
        initial_activations=None,

        enable_myo=settings["enable_myo"],
        enable_shear=settings["enable_shear"],
        enable_meta=settings["enable_meta"],

        dt=DT,
        max_iterations=MAX_ITERATIONS,
    )

    if not control_ar["converged"]:

        raise RuntimeError(
            f"{mechanism_name}: "
            f"failed to converge at "
            f"{CONTROL_PRESSURE:.0f} mmHg."
        )

    rows = []

    rows.append(
        extract_row(
            CONTROL_PRESSURE,
            control_ar,
            mechanism_name,
        )
    )

    print()
    print("100 mmHg equilibrium")

    print_result(
        CONTROL_PRESSURE,
        control_ar,
    )

    # ========================================================
    # DOWNWARD SWEEP
    #
    # 100 -> 90 -> 80 -> ... -> 20
    # ========================================================

    print()
    print("DOWNWARD SWEEP")

    # Start from THIS mechanism's
    # 100-mmHg equilibrium.

    current_data = (
        control_ar["data"].copy()
    )

    current_activations = {

        "Large Arteriole":
            control_ar["activation_LA"],

        "Small Arteriole":
            control_ar["activation_SA"],
    }

    # Myo + Shear is evaluated from 40 to 200 mmHg.
    # The other mechanism states retain the original 20 to 200 mmHg range.
    min_pressure = (
        40.0
        if mechanism_name == "Myo + Shear"
        else MIN_PRESSURE
    )

    downward_pressures = np.arange(
        CONTROL_PRESSURE - PRESSURE_STEP,
        min_pressure - PRESSURE_STEP,
        -PRESSURE_STEP
    )

    for pressure in downward_pressures:

        result = solve_pressure(
            pressure,
            current_data,
            current_activations,
            control_results,
            parameters,
            settings,
        )

        rows.append(
            extract_row(
                pressure,
                result,
                mechanism_name,
            )
        )

        print_result(
            pressure,
            result,
        )

        if not result["converged"]:

            print(
                f"{mechanism_name}: "
                "downward sweep stopped."
            )

            break

        # ----------------------------------------
        # ONLY continuation inside THIS mechanism
        # ----------------------------------------

        current_data = (
            result["data"].copy()
        )

        current_activations = {

            "Large Arteriole":
                result["activation_LA"],

            "Small Arteriole":
                result["activation_SA"],
        }

    # ========================================================
    # UPWARD SWEEP
    #
    # IMPORTANT:
    # Do NOT continue from the 20-mmHg endpoint.
    #
    # Restart from THIS mechanism's
    # original 100-mmHg equilibrium.
    # ========================================================

    print()
    print("UPWARD SWEEP")

    current_data = (
        control_ar["data"].copy()
    )

    current_activations = {

        "Large Arteriole":
            control_ar["activation_LA"],

        "Small Arteriole":
            control_ar["activation_SA"],
    }

    upward_pressures = np.arange(
        CONTROL_PRESSURE + PRESSURE_STEP,
        MAX_PRESSURE + PRESSURE_STEP,
        PRESSURE_STEP
    )

    for pressure in upward_pressures:

        result = solve_pressure(
            pressure,
            current_data,
            current_activations,
            control_results,
            parameters,
            settings,
        )

        rows.append(
            extract_row(
                pressure,
                result,
                mechanism_name,
            )
        )

        print_result(
            pressure,
            result,
        )

        if not result["converged"]:

            print(
                f"{mechanism_name}: "
                "upward sweep stopped."
            )

            break

        # ----------------------------------------
        # ONLY continuation inside THIS mechanism
        # ----------------------------------------

        current_data = (
            result["data"].copy()
        )

        current_activations = {

            "Large Arteriole":
                result["activation_LA"],

            "Small Arteriole":
                result["activation_SA"],
        }

    # ========================================================
    # DATAFRAME
    # ========================================================

    results = pd.DataFrame(rows)

    results = (
        results
        .sort_values("pressure_mmhg")
        .reset_index(drop=True)
    )

    # ========================================================
    # NORMALIZATION
    # ========================================================

    control = results.loc[
        np.isclose(
            results["pressure_mmhg"],
            CONTROL_PRESSURE
        )
    ]

    if control.empty:

        raise RuntimeError(
            f"{mechanism_name}: "
            "100 mmHg result missing."
        )

    Q100 = (
        control.iloc[0]["flow_ul_min"]
    )

    results["normalized_perfusion"] = (
        results["flow_ul_min"]
        / Q100
    )

    # ========================================================
    # SAVE THIS STATE INDIVIDUALLY
    # ========================================================

    safe_name = (
        mechanism_name
        .lower()
        .replace(" + ", "_")
        .replace(" ", "_")
    )

    results.to_csv(
        f"results_{safe_name}.csv",
        index=False
    )

    print()
    print(
        f"Finished mechanism: "
        f"{mechanism_name}"
    )

    print(
        "================================================"
    )

    # Function exits here.
    # All local solver states disappear.
    # Next mechanism starts from INPUT_FILE again.

    return results


# ============================================================
# MAIN
# ============================================================

def main():

    all_results = []

    # ========================================================
    # RUN THREE STATES SEQUENTIALLY
    #
    # Each call performs a complete hard reset.
    # ========================================================

    for (
        mechanism_name,
        settings
    ) in MECHANISM_STATES.items():

        result = run_pressure_sweep(
            mechanism_name,
            settings,
        )

        all_results.append(
            result.copy()
        )

        # There is deliberately no transfer of:
        #
        # data
        # activation
        # control_results
        # parameters
        #
        # into the next mechanism.

    # ========================================================
    # COMBINE ONLY THE FINISHED DATAFRAMES
    # ========================================================

    results_all = pd.concat(
        all_results,
        ignore_index=True
    )

    results_all.to_csv(
        "results_pressure_sweep_three_states.csv",
        index=False
    )

    # ========================================================
    # PRINT RESULTS
    # ========================================================

    print()
    print("================================================")
    print(" FINAL COMBINED RESULTS")
    print("================================================")

    print(
        results_all[
            [
                "mechanism",
                "pressure_mmhg",
                "flow_ul_min",
                "normalized_perfusion",
                "diameter_LA_um",
                "activation_LA",
                "diameter_SA_um",
                "activation_SA",
                "iterations",
                "converged",
            ]
        ].sort_values(
            [
                "mechanism",
                "pressure_mmhg"
            ]
        ).to_string(
            index=False
        )
    )

    # ========================================================
    # THREE STATES ON ONE PLOT
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    plot_styles = {

        "All": {
            "linestyle": "-",
            "marker": "o",
        },

        "Myo + Shear": {
            "linestyle": "--",
            "marker": "s",
        },

        "Myo": {
            "linestyle": "-.",
            "marker": "^",
        },
    }

    for mechanism_name in MECHANISM_STATES:

        subset = results_all.loc[
            results_all["mechanism"]
            == mechanism_name
        ].copy()

        subset = subset.sort_values(
            "pressure_mmhg"
        )

        ax.plot(
            subset["pressure_mmhg"],
            subset["normalized_perfusion"],

            linestyle=
                plot_styles[
                    mechanism_name
                ]["linestyle"],

            marker=
                plot_styles[
                    mechanism_name
                ]["marker"],

            color="black",

            markersize=5,

            label=mechanism_name,
        )

    ax.axhline(
        1.0,
        linestyle=":",
        color="black",
        linewidth=1,
    )

    ax.set_xlabel(
        "Arterial pressure (mmHg)"
    )

    ax.set_ylabel(
        "Normalized perfusion"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        "pressure_vs_normalized_perfusion_three_states.png",
        dpi=300
    )
        # Save selected sweep results
    results_sweep = results_all[
        [
            "mechanism",
            "pressure_mmhg",
            "flow_ul_min",
            "normalized_perfusion",
            "diameter_LA_um",
            "diameter_SA_um",
            "activation_LA",
            "activation_SA",
        ]
    ].copy()

    results_sweep = results_sweep.sort_values(
        [
            "mechanism",
            "pressure_mmhg",
        ]
    ).reset_index(drop=True)

    results_sweep.to_csv(
        "results_sweep.csv",
        index=False,
    )
    
    plt.show()


if __name__ == "__main__":
    main()
