
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import gc

# Configuration
DATA = Path(".")
OUT = Path("plots")

SELECTED_IC = 1

N_FLUX = 5
N_TIME = 2

if SELECTED_IC not in (1, 2, 3):
    raise ValueError("SELECTED_IC must be 1, 2, or 3")

OUTPUT_INTERVAL = 0.0002
TEND = 0.002

EXACT_TIMES = np.arange(
    0.0, TEND + OUTPUT_INTERVAL/2, OUTPUT_INTERVAL
)

VARIABLES = {
    "rho": {
        "column": 1,
        "label": r"$\rho$",
        "title": "Density",
    },
    "u": {
        "column": 2,
        "label": r"$u$",
        "title": "Velocity",
    },
    "p": {
        "column": 3,
        "label": r"$p$",
        "title": "Pressure",
    },
    "T": {
        "column": 4,
        "label": r"$T$",
        "title": "Temperature",
    },
}

FLUX_NAMES = {
    1: "Steger-Warming",
    2: "Van Leer (Local)",
    3: "Van Leer (Common)",
    4: "AUSM",
    5: "Roe",
}

TIME_NAMES = {
    1: "Forward Euler",
    2: "Runge-Kutta 4",
}

# Plot style
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.labelsize": 18,
    "axes.titlesize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 11,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,
    "ytick.right": False,
    "xtick.major.size": 7,
    "ytick.major.size": 7,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
})

# Directory structure
OUT.mkdir(parents=True, exist_ok=True)

# Formatting
def format_axes(ax, ylabel):

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(ylabel)

    ax.grid(False)

    ax.tick_params(
        which="major",
        direction="in",
        length=7,
        width=1,
        pad=7,
    )

    ax.tick_params(
        which="minor",
        direction="in",
        length=4,
        width=0.8,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("black")

# Read blank-line-separated solution blocks
def read_blocks(filename, ncols=5):

    blocks = []
    block = []

    with open(filename, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                if block:
                    blocks.append(np.asarray(block, dtype=float))
                    block = []
                continue

            if line.startswith("#"):
                continue

            values = line.split()

            if len(values) != ncols:
                continue

            block.append([float(v) for v in values])

    if block:
        blocks.append(np.asarray(block, dtype=float))

    if not blocks:
        raise ValueError(f"No data blocks found in {filename}")

    return blocks

# Read numerical and exact solutions
def load_case(ic, flux, time_integrator):

    numerical_file = (DATA /f"results-IC{ic}-F{flux}-T{time_integrator}.dat")
    exact_file = DATA / f"exact-IC{ic}.dat"

    if not numerical_file.exists():
        raise FileNotFoundError(numerical_file)

    if not exact_file.exists():
        raise FileNotFoundError(exact_file)

    numerical = read_blocks(numerical_file)
    exact = read_blocks(exact_file)

    numerical_times = np.arange(
        len(numerical)
    ) * OUTPUT_INTERVAL

    if len(exact) != len(EXACT_TIMES):
        raise ValueError(
            f"{exact_file}: expected {len(EXACT_TIMES)} "
            f"exact blocks, found {len(exact)}"
        )

    if len(numerical) != len(numerical_times):
        raise ValueError("Unexpected numerical block count")

    return numerical, exact, numerical_times

# Plot numerical vs exact solution
def plot_solution(ic, flux, time_integrator,time, numerical, exact):

    # Output directory
    directory = (OUT / f"IC{ic}" / f"F{flux}" / f"T{time_integrator}")
    directory.mkdir(parents=True, exist_ok=True)

    # Create 2x2 figure
    fig, axes = plt.subplots(
        2, 2,
        figsize=(11, 8),
    )

    axes = axes.flatten()

    # Subfigure labels in report order
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for ax, (variable, info), panel in zip(axes, VARIABLES.items(), panel_labels):

        col = info["column"]

        # Numerical solution
        ax.plot(
            numerical[:, 0],
            numerical[:, col],
            color="blue",
            linewidth=1.8,
            label="Numerical",
        )

        # Exact solution
        ax.plot(
            exact[:, 0],
            exact[:, col],
            color="black",
            linestyle="--",
            linewidth=1.6,
            label="Exact",
        )

        format_axes(ax, info["label"])
        ax.set_title("")

        ax.legend(
            loc="best",
            frameon=True,
            fancybox=False,
            framealpha=1,
            edgecolor="black",
        )

        # Panel label below each subplot
        ax.text(
            0.5,
            -0.28,
            panel,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=15,
        )

    fig.tight_layout()

    # Leave sufficient room for panel labels
    fig.subplots_adjust(
        hspace=0.55,
        wspace=0.30,
        bottom=0.12,
    )

    filename = directory / f"t{time:.4f}.png"

    fig.savefig(
        filename,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

    del fig, axes

# def plot_solution(ic, flux, time_integrator,
#                   time, numerical, exact):

#     # Output directory
#     directory = (
#         OUT /
#         f"IC{ic}" /
#         f"F{flux}" /
#         f"T{time_integrator}"
#     )

#     directory.mkdir(parents=True, exist_ok=True)

#     # Create 2x2 figure
#     fig, axes = plt.subplots(
#         2, 2,
#         figsize=(11, 8),
#     )

#     axes = axes.flatten()

#     for ax, (variable, info) in zip(
#         axes, VARIABLES.items()
#     ):

#         col = info["column"]

#         # Numerical solution
#         ax.plot(
#             numerical[:, 0],
#             numerical[:, col],
#             color="blue",
#             linewidth=1.8,
#             label="Numerical",
#         )

#         # Exact solution
#         ax.plot(
#             exact[:, 0],
#             exact[:, col],
#             color="black",
#             linestyle="--",
#             linewidth=1.6,
#             label="Exact",
#         )

#         format_axes(ax, info["label"])

#         ax.set_title(info["title"])

#         ax.legend(
#             loc="best",
#             frameon=True,
#             fancybox=False,
#             framealpha=1,
#             edgecolor="black",
#         )

#     fig.suptitle(
#         rf"IC {ic} | {FLUX_NAMES[flux]} | "
#         rf"{TIME_NAMES[time_integrator]} | "
#         rf"$t={time:.4f}$",
#         fontsize=15,
#         y=0.99,
#     )

#     fig.tight_layout(rect=[0, 0, 1, 0.96])

#     filename = directory / f"t{time:.4f}.png"

#     fig.savefig(
#         filename,
#         dpi=600,
#         bbox_inches="tight",
#     )

#     plt.close(fig)

#     # Explicitly release references
#     del fig, axes, ax


# Main post-processing loop
# for ic in range(1, N_IC + 1):

#     print(f"\nProcessing IC {ic}")

#     for flux in range(1, N_FLUX + 1):

#         for time_integrator in range(1, N_TIME + 1):

#             print(
#                 f"  Flux: {flux} ({FLUX_NAMES[flux]}) | "
#                 f"Time integrator: {time_integrator} "
#                 f"({TIME_NAMES[time_integrator]})"
#             )

#             try:
#                 numerical, exact, times = load_case(
#                     ic,
#                     flux,
#                     time_integrator,
#                 )

#                 for i, (time, data) in enumerate(
#                     zip(times, numerical)
#                 ):

#                     # Find corresponding exact snapshot
#                     exact_index = np.argmin(
#                         np.abs(
#                             np.asarray(EXACT_TIMES) - time
#                         )
#                     )

#                     if not np.isclose(
#                         EXACT_TIMES[exact_index],
#                         time,
#                         atol=1e-9,
#                     ):
#                         print(
#                             f"    Skipping t={time:.4f}: "
#                             "no matching exact snapshot"
#                         )
#                         continue

#                     plot_solution(
#                         ic,
#                         flux,
#                         time_integrator,
#                         time,
#                         data,
#                         exact[exact_index],
#                     )

#                     gc.collect()

#                 del numerical, exact, times
#                 gc.collect()

#             except (FileNotFoundError, ValueError) as error:

#                 print(f"    WARNING: {error}")
#                 continue

# print("\nPost-processing complete.")

# Main post-processing loop
ic = SELECTED_IC

print(f"\nProcessing IC {ic}")

for flux in range(1, N_FLUX + 1):

    for time_integrator in range(1, N_TIME + 1):

        print(
            f"  Flux: {flux} ({FLUX_NAMES[flux]}) | "
            f"Time integrator: {time_integrator} "
            f"({TIME_NAMES[time_integrator]})"
        )

        try:
            numerical, exact, times = load_case(
                ic,
                flux,
                time_integrator,
            )

            for time, data in zip(times, numerical):

                # Find corresponding exact snapshot
                exact_index = np.argmin(
                    np.abs(np.asarray(EXACT_TIMES) - time)
                )

                if not np.isclose(
                    EXACT_TIMES[exact_index],
                    time,
                    atol=1e-9,
                ):
                    print(
                        f"    Skipping t={time:.4f}: "
                        "no matching exact snapshot"
                    )
                    continue

                plot_solution(
                    ic,
                    flux,
                    time_integrator,
                    time,
                    data,
                    exact[exact_index],
                )

            del numerical, exact, times
            gc.collect()

        except (FileNotFoundError, ValueError) as error:

            print(f"    WARNING: {error}")
            continue

print("\nPost-processing complete.")