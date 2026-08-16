import math
import pandas as pd

from inputs import INPUT_FILE


HD = 0.45


def c_value(D):
    return (
        (0.8 + math.exp(-0.075 * D))
        * (-1 + 1 / (1 + 1e-11 * D**12))
        + 1 / (1 + 1e-11 * D**12)
    )


def viscosity_045(D):
    return (
        6 * math.exp(-0.085 * D)
        + 3.2
        - 2.44 * math.exp(-0.06 * D**0.645)
    )


def viscosity(D):

    C = c_value(D)
    mu_045 = viscosity_045(D)

    hematocrit = (
        ((1 - HD)**C - 1)
        / ((1 - 0.45)**C - 1)
    )

    correction = (
        D / (D - 1.1)
    )**2

    return (
        1
        + (mu_045 - 1)
        * hematocrit
        * correction
    ) * correction


def main():

    data = pd.read_csv(INPUT_FILE)

    results = data[
        ["compartment", "diameter_um"]
    ].copy()

    results["viscosity_cp"] = (
        results["diameter_um"].apply(viscosity)
    )

    results.to_csv(
        "viscosity_results.csv",
        index=False
    )

    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
