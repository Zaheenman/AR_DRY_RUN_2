from pathlib import Path

import matplotlib.pyplot as plt


PLOT_DIR = Path(__file__).resolve().parent / "plots"
PLOT_DIR.mkdir(exist_ok=True)


def add_boundaries(profile):

    names = [
        "A",
        "LA",
        "SA",
        "C",
        "SV",
        "LV",
    ]

    boundaries = profile["boundaries"]

    for x in boundaries:
        plt.axvline(
            x * 100,
            linewidth=0.5,
            alpha=0.3
        )

    y_min, _ = plt.ylim()

    for i, name in enumerate(names):

        x1 = boundaries[i] * 100
        x2 = boundaries[i + 1] * 100

        plt.text(
            (x1 + x2) / 2,
            y_min,
            name,
            ha="center",
            va="bottom"
        )


def plot_figure2(profiles):

    styles = {
        1.0: ":",
        8.28: "--",
        20.0: "-",
    }

    labels = {
        1.0: "M0 = 1",
        8.28: "M0 = 8.28",
        20.0: "M0 = 20",
    }

    plot_saturation(
        profiles,
        styles,
        labels
    )

    plot_atp(
        profiles,
        styles,
        labels
    )

    plot_scr(
        profiles,
        styles,
        labels
    )


def plot_saturation(
    profiles,
    styles,
    labels
):

    plt.figure()

    for M0, profile in profiles.items():

        plt.plot(
            profile["x"] * 100,
            profile["S"],
            styles[M0],
            label=labels[M0]
        )

    plt.xlabel(
        "Distance along vascular pathway (cm)"
    )

    plt.ylabel(
        "Oxygen saturation"
    )

    plt.ylim(
        0,
        1
    )

    add_boundaries(
        profiles[8.28]
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / "figure2_saturation.png",
        dpi=300
    )

    plt.close()


def plot_atp(
    profiles,
    styles,
    labels
):

    plt.figure()

    for M0, profile in profiles.items():

        plt.plot(
            profile["x"] * 100,
            profile["C"] * 1000,
            styles[M0],
            label=labels[M0]
        )

    plt.xlabel(
        "Distance along vascular pathway (cm)"
    )

    plt.ylabel(
        "ATP concentration (µM)"
    )

    add_boundaries(
        profiles[8.28]
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / "figure2_ATP.png",
        dpi=300
    )

    plt.close()


def plot_scr(
    profiles,
    styles,
    labels
):

    plt.figure()

    for M0, profile in profiles.items():

        plt.plot(
            profile["x"] * 100,
            profile["Scr"] * 1e5,
            styles[M0],
            label=labels[M0]
        )

    reference = profiles[8.28]

    x_LA = reference["x_LA"] * 100
    x_SA = reference["x_SA"] * 100

    plt.axvline(
        x_LA,
        linestyle=":",
        linewidth=1
    )

    plt.axvline(
        x_SA,
        linestyle=":",
        linewidth=1
    )

    y_max = plt.ylim()[1]

    plt.text(
        x_LA,
        y_max * 0.95,
        "LA",
        ha="center"
    )

    plt.text(
        x_SA,
        y_max * 0.95,
        "SA",
        ha="center"
    )

    plt.xlabel(
        "Distance along vascular pathway (cm)"
    )

    plt.ylabel(
        "Scr (µM·cm)"
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR / "figure2_Scr.png",
        dpi=300
    )

    plt.close()


def plot_figure4(results):

    plot_activation_diameter(
        results,
        "LA",
        "Large Arteriole"
    )

    plot_activation_diameter(
        results,
        "SA",
        "Small Arteriole"
    )

    plot_tone_components(
        results,
        "LA",
        "Large Arteriole"
    )

    plot_tone_components(
        results,
        "SA",
        "Small Arteriole"
    )

    plot_ar_flow(
        results
    )


def plot_activation_diameter(
    results,
    vessel,
    title
):

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()

    activation = ax1.plot(
        results["M0"],
        results[f"activation_{vessel}"],
        "-",
        label="Activation"
    )

    diameter = ax2.plot(
        results["M0"],
        results[f"diameter_{vessel}_um"],
        "--",
        label="Diameter"
    )

    ax1.set_xlabel(
        r"Consumption "
        r"(cm$^3$ O$_2$/100 cm$^3$/min)"
    )

    ax1.set_ylabel(
        "Activation"
    )

    ax2.set_ylabel(
        "Diameter (µm)"
    )

    ax1.set_xlim(
        0,
        25
    )

    ax1.set_ylim(
        0,
        1
    )

    ax1.set_title(
        title
    )

    lines = (
        activation
        + diameter
    )

    ax1.legend(
        lines,
        [
            line.get_label()
            for line in lines
        ],
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        PLOT_DIR
        / f"figure4_{vessel}_activation_diameter.png",
        dpi=300
    )

    plt.close(
        fig
    )


def plot_tone_components(
    results,
    vessel,
    title
):

    myogenic = results[
        f"{vessel}_myogenic"
    ]

    shear = abs(
        results[
            f"{vessel}_shear"
        ]
    )

    metabolic = abs(
        results[
            f"{vessel}_metabolic"
        ]
    )

    plt.figure()

    plt.plot(
        results["M0"],
        myogenic,
        "-",
        label=r"$C_{myo}T$ (+)"
    )

    plt.plot(
        results["M0"],
        shear,
        "--",
        label=r"$C_{shear}\tau$ (-)"
    )

    plt.plot(
        results["M0"],
        metabolic,
        "-.",
        label=r"$C_{meta}S_{CR}$ (-)"
    )

    plt.xlabel(
        r"Consumption "
        r"(cm$^3$ O$_2$/100 cm$^3$/min)"
    )

    plt.ylabel(
        r"Components of $S_{tone}$"
    )

    plt.title(
        title
    )

    plt.xlim(
        0,
        25
    )

    plt.ylim(
        0,
        16
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        PLOT_DIR
        / f"figure4_{vessel}_tone_components.png",
        dpi=300
    )

    plt.close()


def plot_ar_flow(results):

    plt.figure()

    plt.plot(
        results["M0"],
        results["flow_ul_min"],
        "-"
    )

    plt.xlabel(
        r"Consumption "
        r"(cm$^3$ O$_2$/100 cm$^3$/min)"
    )

    plt.ylabel(
        "Total flow (µL/min)"
    )

    plt.xlim(
        0,
        25
    )

    plt.tight_layout()

    plt.savefig(
        PLOT_DIR
        / "figure4_total_flow.png",
        dpi=300
    )

    plt.close()
