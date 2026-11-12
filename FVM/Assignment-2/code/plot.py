from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Style
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.labelsize": 18,
    "axes.titlesize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
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


# Paths
DATA = Path(".")
OUT = Path("plots")

for d in [
    OUT/"ic",
    OUT/"sol"/"t01",
    OUT/"sol"/"t05",
    OUT/"res"/"t01",
    OUT/"res"/"t05",
]:
    d.mkdir(parents=True, exist_ok=True)

# Common plotting function
def format_axes(ax, ylabel):
    ax.set_xlabel(r"$x$", labelpad=10)
    ax.set_ylabel(ylabel, labelpad=10)
    ax.grid(False)

    ax.tick_params(which="major", direction="in", length=7, width=1.0, pad=7)
    ax.tick_params(which="minor", direction="in", length=4, width=0.8)

    for s in ax.spines.values():
        s.set_linewidth(1.2)
        s.set_color("black")

def finish(fig, ax, filename, legend_specific=False):
    if legend_specific:
        ax.legend(
            loc="center left",
            bbox_to_anchor=(0.65, 0.3),
            frameon=True,
            fancybox=False,
            framealpha=1,
            edgecolor="black",
        )
        fig.tight_layout(rect=[0, 0, 0.84, 1])

    else:
        leg = ax.legend(
            loc="best",
            frameon=True,
            fancybox=False,
            framealpha=1,
            edgecolor="black",
        )

        fig.canvas.draw()

        bbox = leg.get_window_extent(fig.canvas.get_renderer())
        bbox = bbox.transformed(ax.transAxes.inverted())

        # Adjust the Legend
        if bbox.x0 < 0.5 and bbox.y0 > 0.5:
            bbox.x0 += 0.02
            bbox.x1 += 0.02
            bbox.y0 -= 0.01
            bbox.y1 -= 0.01
            leg.set_bbox_to_anchor(bbox, transform=ax.transAxes)

        if bbox.x0 < 0.5 and bbox.y0 < 0.5:
            bbox.x0 += 0.02
            bbox.x1 += 0.02
            bbox.y0 += 0.05
            bbox.y1 += 0.05
            leg.set_bbox_to_anchor(bbox, transform=ax.transAxes)

        fig.tight_layout()
        
    fig.savefig(filename, dpi=600, bbox_inches="tight")
    plt.close(fig)

# Read results
def read_results(filename):

    blocks = []
    block = []

    with open(filename) as f:
        for line in f:

            line = line.strip()

            if not line:
                if block:
                    blocks.append(np.array(block))
                    block = []
                continue

            if line.startswith("#"):
                continue

            vals = line.split()

            if len(vals) == 7:
                block.append([float(v) for v in vals])

        if block:
            blocks.append(np.array(block))

    return blocks[0], blocks[1]

# Initial conditions
def initial_conditions():
    x = np.linspace(-1, 1, 21)
    xc = 0.5 * (x[1:] + x[:-1])
    u = {
        1: np.where(xc <= 0, 1, 0),
        2: np.where(xc <= 0, 0, 1),
        3: np.where(xc <= 0, -1, 1),
        4: np.where(xc <= 0, 1, -1),
        5: np.where(
            xc <= -0.5,
            1,
            np.where(xc >= 0.5, 0, 0.5 - xc),
        ),
    }
    return xc, u

def plot_initial_conditions():

    x, u = initial_conditions()
    colors = ["red", "blue", "green", "black", "magenta"]

    fig, ax = plt.subplots(figsize=(7, 5))

    for ic in range(1, 6):
        ax.plot(
            x, u[ic],
            color=colors[ic - 1],
            linewidth=1.8,
            marker="o",
            markersize=3.5,
            label=rf"Case {ic}",
        )

    format_axes(ax, r"$u(x, t = 0)$")

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xticks(np.arange(-1, 1.1, 0.5))
    ax.set_yticks(np.arange(-1, 1.1, 0.5))

    finish(fig, ax, OUT/"ic"/"ic.png", legend_specific=True)


# Exact analytical solutions for Burgers' equation
def exact_solution(ic, x, t):

    x = np.asarray(x)

    # Case 1: u_L = 1, u_R = 0
    # Moving shock with speed s = 1/2
    
    if ic == 1:
        return np.where(x < 0.5 * t, 1.0, 0.0)

    # Case 2: u_L = 0, u_R = 1
    # Rarefaction
    
    elif ic == 2:
        u = np.empty_like(x)

        u[x < 0] = 0.0
        u[x > t] = 1.0

        mask = (x >= 0) & (x <= t)
        u[mask] = x[mask] / t

        return u

    # Case 3: u_L = -1, u_R = 1
    # Rarefaction
    
    elif ic == 3:
        u = np.empty_like(x)

        u[x < -t] = -1.0
        u[x > t] = 1.0

        mask = (x >= -t) & (x <= t)
        u[mask] = x[mask] / t

        return u

    # Case 4: u_L = 1, u_R = -1
    # Stationary shock
    
    elif ic == 4:
        return np.where(x < 0.0, 1.0, -1.0)

    # Case 5:
    #
    # u = 1             x <= -0.5
    # u = 0.5 - x       -0.5 < x < 0.5
    # u = 0             x >= 0.5
    #
    # Valid for t < 1, before shock formation.
    
    elif ic == 5:
        u = np.empty_like(x)

        x_left = -0.5 + t
        x_right = 0.5

        u[x < x_left] = 1.0
        u[x > x_right] = 0.0

        mask = (x >= x_left) & (x <= x_right)
        u[mask] = (0.5 - x[mask]) / (1.0 - t)

        return u

    else:
        raise ValueError(f"Unknown IC case: {ic}")


# Numerical + exact solution / residual plotting
def plot_result(ic, time, data, residual=False):
    x = data[:, 0]
    if residual:
        columns = [2, 4, 6]
        ylabel = r"$R$"
        directory = OUT/"res"/f"t{int(time * 10):02d}"
    else:
        columns = [1, 3, 5]
        ylabel = r"$u$"
        directory = OUT/"sol"/f"t{int(time * 10):02d}"
    fig, ax = plt.subplots(figsize=(7, 5))

    # FVS
    ax.plot(
        x, data[:, columns[0]],
        color="blue",
        linestyle="-",
        linewidth=2,
        label="FVS",
    )

    # FDS
    ax.plot(
        x, data[:, columns[1]],
        color="red",
        linestyle="--",
        linewidth=2,
        label="FDS",
    )
    
    # Godunov
    ax.plot(
        x, data[:, columns[2]],
        color="black",
        linestyle="None",
        marker="o",
        markersize=4,
        label="Godunov",
    )

    # Analytical exact solution
    if not residual:

        x_exact = np.linspace(x.min(),x.max(),21)
        u_exact = exact_solution(ic,x_exact,time)

        ax.plot(
            x_exact,
            u_exact,
            color="black",
            linestyle="None",
            linewidth=1.8,
            marker="s",
            markersize=4.5,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=1.2,
            label="Exact",
            zorder=2,
        )

    format_axes(ax, ylabel)
    finish(fig,ax,directory/f"IC{ic}.png")

# Main
plot_initial_conditions()

for ic in range(1, 6):
    filename = DATA/f"results-{ic}.dat"

    if not filename.exists():
        print(f"Skipping IC {ic}: {filename} not found")
        continue

    print(f"Processing {filename}")

    t01, t05 = read_results(filename)

    for time, data in [(0.1, t01), (0.5, t05)]:
        plot_result(ic, time, data, residual=False)
        plot_result(ic, time, data, residual=True)

print("\nDone.")