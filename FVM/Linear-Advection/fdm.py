#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
c = 1.0
N = 100
x = np.arange(N)/N
dx = 1.0/N
T = 1.0

mus = [1.0]    # [0.25, 0.5, 0.75, 1.0, 1.25]
m = 1

# Initial conditions
def initial_condition_a(x):
    return np.sin(2*np.pi*m*x)

def initial_condition_b(x):
    u = np.zeros_like(x)
    u[(x >= 0.4) & (x <= 0.6)] = 1.0
    return u

# Exact solution
def exact_solution(x, t, ic):
    if ic == "sine":
        return np.sin(2*np.pi*m*(x - c*t))
    shift = int(round(c*t/dx))
    return np.roll(initial_condition_b(x), shift)

# Numerical schemes
def FTBS(u, mu):
    return u - mu*(u - np.roll(u, 1))

def LaxWendroff(u, mu):
    return (u - 0.5*mu*(np.roll(u, -1) - np.roll(u, 1)) + 0.5*mu**2*(np.roll(u, -1) - 2*u + np.roll(u, 1)))

def FTCS(u, mu):
    return u - 0.5*mu*(np.roll(u, -1) - np.roll(u, 1))

schemes = {"bs": FTBS, "lw": LaxWendroff, "cs": FTCS}

# Solver
def solve(u0, scheme, mu):
    dt = mu*dx/c
    nt = int(T/dt)
    u = u0.copy()
    data = [u.copy()]

    for _ in range(nt):
        u = scheme(u, mu)
        data.append(u.copy())

    return np.array(data), dt

# Animation
def animate(u0, scheme, mu, ic):
    data, dt = solve(u0, scheme, mu)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.grid()

    numerical, = ax.plot(x, data[0], lw=2, label="Numerical")

    exact, = ax.plot(x, exact_solution(x, 0, ic), "k--", label="Exact")
    ax.legend()

    def update(n):
        t = n*dt
        numerical.set_ydata(data[n])
        exact.set_ydata(exact_solution(x, t, ic))
        ax.set_title(f"{scheme.__name__} | μ={mu}")
        return numerical, exact

    ani = FuncAnimation(fig, update, frames=len(data), interval=30, blit=True, repeat=True)

    plt.tight_layout()
    plt.show()

# Main execution
if __name__ == "__main__":

    print("\nLinear Advection Equation \n")

    ic = input("Initial condition [sine/square]: ").strip().lower()
    scheme = input("Scheme [bs/lw/cs]: ").strip().lower()

    if ic not in ("sine", "square"):
        raise ValueError("Initial condition must be 'sine' or 'square'.")

    if scheme not in schemes:
        raise ValueError("Scheme must be 'bs', 'lw', or 'cs'.")

    u0 = initial_condition_a(x) if ic == "sine" else initial_condition_b(x)

    for mu in mus:
        animate(u0, schemes[scheme], mu, ic)
