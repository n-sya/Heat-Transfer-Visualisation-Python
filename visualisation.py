import os

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from variables import wire_radius

# Output directory
OUTPUT_DIRECTORY = "outputs"


# Create output directory
def create_output_directory():
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)


# Convert polar coordinates to Cartesian coordinates
def polar_to_cartesian(radius, theta):
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    return x, y


# Plot radial temperature distribution
def plot_temperature_profile(T, r, solver_name):
    create_output_directory()

    r_full = np.concatenate((-r[::-1][:-1], r))
    T_full = np.concatenate((T[-1, ::-1][:-1], T[-1, :]))

    fig, ax = plt.subplots()
    ax.plot(r_full * 1000, T_full)
    ax.grid()
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"{solver_name} Temperature Distribution")

    fig.tight_layout()
    fig.savefig(
        os.path.join(
            OUTPUT_DIRECTORY,
            f"{solver_name.lower()}_temperature_profile.png",
        ),
        dpi=300,
    )

    plt.show()


# Plot final temperature heatmap
def plot_temperature_heatmap(T, r, solver_name):
    create_output_directory()

    theta = np.linspace(0, 2 * np.pi, 360)
    theta_edges = np.linspace(0, 2 * np.pi, len(theta) + 1)
    radius_edges = np.linspace(0, wire_radius, len(r) + 1)

    radius_grid, theta_grid = np.meshgrid(radius_edges, theta_edges)
    x, y = polar_to_cartesian(radius_grid, theta_grid)

    T_final = T[-1, :]
    T_assign = np.tile(T_final, (len(theta_edges) - 1, 1))

    fig, ax = plt.subplots()

    heatmap = ax.pcolormesh(
        x * 1000,
        y * 1000,
        T_assign,
        shading="auto",
    )

    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"{solver_name} Temperature Distribution")
    ax.axis("equal")

    fig.colorbar(heatmap, ax=ax, label="Temperature (K)")
    fig.tight_layout()

    fig.savefig(
        os.path.join(
            OUTPUT_DIRECTORY,
            f"{solver_name.lower()}_temperature_heatmap.png",
        ),
        dpi=300,
    )

    plt.show()


# Animate temperature distribution
def animate_temperature(T, r, t, solver_name, frame_interval=400):
    create_output_directory()

    theta = np.linspace(0, 2 * np.pi, 360)
    theta_edges = np.linspace(0, 2 * np.pi, len(theta) + 1)
    radius_edges = np.linspace(0, wire_radius, len(r) + 1)

    radius_grid, theta_grid = np.meshgrid(radius_edges, theta_edges)
    x, y = polar_to_cartesian(radius_grid, theta_grid)

    fig, ax = plt.subplots()

    initial_temperature = np.tile(
        T[0],
        (len(theta_edges) - 1, 1),
    )

    heatmap = ax.pcolormesh(
        x * 1000,
        y * 1000,
        initial_temperature,
        shading="auto",
        vmin=np.min(T),
        vmax=np.max(T),
    )

    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.axis("equal")

    fig.colorbar(heatmap, ax=ax, label="Temperature (K)")

    frames = range(0, len(t), frame_interval)

    def update(frame):
        temperature = np.tile(
            T[frame],
            (len(theta_edges) - 1, 1),
        )

        heatmap.set_array(temperature.ravel())
        ax.set_title(f"{solver_name} Temperature at t = {t[frame]:.2f} s")

        return [heatmap]

    temperature_animation = animation.FuncAnimation(
        fig,
        update,
        frames=frames,
        blit=False,
        repeat=False,
    )

    temperature_animation.save(
        os.path.join(
            OUTPUT_DIRECTORY,
            f"{solver_name.lower()}_temperature_animation.mp4",
        ),
        writer="ffmpeg",
    )

    plt.close(fig)


# Compare final temperature profiles
def plot_solver_comparison(T_explicit, T_crank_nicolson, r):
    create_output_directory()

    r_full = np.concatenate((-r[::-1][:-1], r))

    explicit_full = np.concatenate(
        (T_explicit[-1, ::-1][:-1], T_explicit[-1, :])
    )
    crank_nicolson_full = np.concatenate(
        (T_crank_nicolson[-1, ::-1][:-1], T_crank_nicolson[-1, :])
    )

    fig, ax = plt.subplots()

    ax.plot(
        r_full * 1000,
        explicit_full,
        label="Explicit",
    )
    ax.plot(
        r_full * 1000,
        crank_nicolson_full,
        label="Crank-Nicolson",
    )

    ax.grid()
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Solver Comparison")
    ax.legend()

    fig.tight_layout()
    fig.savefig(
        os.path.join(
            OUTPUT_DIRECTORY,
            "solver_comparison.png",
        ),
        dpi=300,
    )

    plt.show()