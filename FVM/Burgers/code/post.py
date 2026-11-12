from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Plot style
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.labelsize": 18, "axes.titlesize": 16,
    "xtick.labelsize": 13, "ytick.labelsize": 13,
    "legend.fontsize": 12, "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": False, "ytick.right": False,
    "xtick.major.size": 7, "ytick.major.size": 7,
    "xtick.minor.size": 4, "ytick.minor.size": 4 })

DATA, OUT = Path("."), Path("plots")

for d in ["ic", "sol/t01", "sol/t05", "res/t01", "res/t05"]:
    (OUT / d).mkdir(parents=True, exist_ok=True)

def format_axes(ax, ylabel):
    ax.set(xlabel=r"$x$", ylabel=ylabel)
    ax.grid(False)
    ax.tick_params(which="major", direction="in", length=7, width=1, pad=7)
    ax.tick_params(which="minor", direction="in", length=4, width=0.8)
    for s in ax.spines.values():
        s.set_linewidth(1.2)
        s.set_color("black")

def finish(fig, ax, filename, legend_specific=False):
    if legend_specific:
        ax.legend(
            loc="center left", bbox_to_anchor=(0.65, 0.3),
            frameon=True, fancybox=False, framealpha=1,
            edgecolor="black",
        )
        fig.tight_layout(rect=[0, 0, 0.84, 1])
    else:
        leg = ax.legend(
            loc="best", frameon=True, fancybox=False,
            framealpha=1, edgecolor="black",
        )
        fig.canvas.draw()
        bbox = leg.get_window_extent(fig.canvas.get_renderer())
        bbox = bbox.transformed(ax.transAxes.inverted())

        if bbox.x0 < 0.5:
            dy = -0.01 if bbox.y0 > 0.5 else 0.05
            bbox.x0 += 0.02
            bbox.x1 += 0.02
            bbox.y0 += dy
            bbox.y1 += dy
            leg.set_bbox_to_anchor(bbox, transform=ax.transAxes)

        fig.tight_layout()

    fig.savefig(filename, dpi=600, bbox_inches="tight")
    plt.close(fig)

def read_results(filename):
    blocks, block = [], []

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line:
                if block:
                    blocks.append(np.array(block))
                    block = []
                continue

            if not line.startswith("#") and len(v := line.split()) == 7:
                block.append([float(x) for x in v])

    if block:
        blocks.append(np.array(block))

    return blocks[0], blocks[1]

def initial_conditions():
    x = np.linspace(-1, 1, 21)
    xc = (x[1:] + x[:-1]) / 2

    u = {
        1: np.where(xc <= 0, 1, 0),
        2: np.where(xc <= 0, 0, 1),
        3: np.where(xc <= 0, -1, 1),
        4: np.where(xc <= 0, 1, -1),
        5: np.where(xc <= -0.5, 1, np.where(xc >= 0.5, 0, 0.5 - xc))}

    return xc, u

def plot_initial_conditions():
    x, u = initial_conditions()
    colors = ["red", "blue", "green", "black", "magenta"]

    fig, ax = plt.subplots(figsize=(7, 5))

    for ic, color in enumerate(colors, 1):
        ax.plot(
            x, u[ic], color=color, linewidth=1.8,
            marker="o", markersize=3.5, label=f"Case {ic}",
        )

    format_axes(ax, r"$u(x, t = 0)$")
    ax.set(
        xlim=(-1, 1), ylim=(-1.2, 1.2),
        xticks=np.arange(-1, 1.1, 0.5),
        yticks=np.arange(-1, 1.1, 0.5),
    )

    finish(fig, ax, OUT / "ic/ic.png", True)

def exact_solution(ic, x, t):
    x = np.asarray(x)

    if ic == 1:
        return np.where(x < t / 2, 1.0, 0.0)

    if ic == 2:
        return np.where(x < 0, 0,
               np.where(x > t, 1, x / t))

    if ic == 3:
        return np.where(x < -t, -1,
               np.where(x > t, 1, x / t))

    if ic == 4:
        return np.where(x < 0, 1.0, -1.0)

    if ic == 5:
        xl = -0.5 + t
        return np.where(x < xl, 1, np.where(x > 0.5, 0, (0.5 - x) / (1 - t)))


def plot_result(ic, time, data, residual=False):
    x = data[:, 0]
    cols = [2, 4, 6] if residual else [1, 3, 5]
    ylabel = r"$R$" if residual else r"$u$"
    directory = OUT / ("res" if residual else "sol") / f"t{int(10*time):02d}"

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(x, data[:, cols[0]], "b-", linewidth=2, label="FVS")
    ax.plot(x, data[:, cols[1]], "r--", linewidth=2, label="FDS")
    ax.plot(x, data[:, cols[2]], "ko", markersize=4, label="Godunov")

    if not residual:
        ax.plot(
            x, exact_solution(ic, x, time), "ks",
            markersize=4.5,
            markerfacecolor="white",
            markeredgewidth=1.2,
            label="Exact",
            zorder=2,
        )

    format_axes(ax, ylabel)
    finish(fig, ax, directory / f"IC{ic}.png")

plot_initial_conditions()

for ic in range(1, 6):
    filename = DATA / f"results-{ic}.dat"

    print(f"Processing {filename}")
    t01, t05 = read_results(filename)

    for time, data in [(0.1, t01), (0.5, t05)]:
        plot_result(ic, time, data)
        plot_result(ic, time, data, True)

print("\nDone.")