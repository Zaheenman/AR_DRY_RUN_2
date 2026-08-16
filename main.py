import numpy as np
import pandas as pd

from inputs import INPUT_FILE
from coupling_solver import solve_figure4
from plots import plot_figure4


def main():

    data = pd.read_csv(
        INPUT_FILE
    )

    M0_values = np.unique(
        np.append(
            np.arange(1.0, 26.0, 1.0),
            8.28
        )
    )

    results = solve_figure4(
        data,
        M0_values
    )

    results.to_csv(
        "results_figure4.csv",
        index=False
    )

    print(
        results[
            [
                "M0",
                "diameter_LA_um",
                "activation_LA",
                "diameter_SA_um",
                "activation_SA",
                "flow_ul_min",
                "iterations",
                "converged",
            ]
        ].to_string(index=False)
    )

    control = results.loc[
        np.isclose(
            results["M0"],
            8.28
        )
    ].iloc[0]

    print()
    print("Control-state check at M0 = 8.28")
    print(
        f"LA: D = {control['diameter_LA_um']:.3f} um, "
        f"A = {control['activation_LA']:.4f}"
    )
    print(
        f"SA: D = {control['diameter_SA_um']:.3f} um, "
        f"A = {control['activation_SA']:.4f}"
    )

    plot_figure4(
        results
    )


if __name__ == "__main__":
    main()
