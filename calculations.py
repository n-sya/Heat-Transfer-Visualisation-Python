import numpy as np

from variables import (
    ambient_temperature,
    convection_coefficient,
    heat_generation,
    initial_temperature,
    rho,
    specific_heat,
    thermal_conductivity,
    wire_radius,
)


# Create radial and temporal grids
def create_grid(n_nodes, dt, t_end, radius=wire_radius):
    dr = radius / (n_nodes - 1)
    r = np.linspace(0, radius, n_nodes)
    t = np.arange(0, t_end + dt, dt)

    return r, t, dr


# Explicit finite difference solver
def explicit_solver(
    n_nodes,
    dt,
    t_end,
    radius=wire_radius,
    density=rho,
    heat_capacity=specific_heat,
    conductivity=thermal_conductivity,
    volumetric_heat_generation=heat_generation,
    initial_temp=initial_temperature,
    ambient_temp=ambient_temperature,
    convection_coefficient_value=convection_coefficient,
):
    thermal_diffusivity = conductivity / (density * heat_capacity)
    source_term = volumetric_heat_generation / (density * heat_capacity)

    r, t, dr = create_grid(n_nodes, dt, t_end, radius)
    n_time_steps = len(t)

    T = np.zeros((n_time_steps, n_nodes))
    T[0, :] = initial_temp

    for n in range(1, n_time_steps):
        previous_T = T[n - 1]

        # Centreline symmetry boundary condition
        T[n, 0] = previous_T[0] + dt * (
            4
            * thermal_diffusivity
            * (previous_T[1] - previous_T[0])
            / dr**2
            + source_term
        )

        # Interior nodes
        for j in range(1, n_nodes - 1):
            radial_conduction = (
                previous_T[j + 1] - 2 * previous_T[j] + previous_T[j - 1]
            ) / dr**2

            cylindrical_term = (previous_T[j + 1] - previous_T[j - 1]) / (
                2 * r[j] * dr
            )

            T[n, j] = previous_T[j] + dt * (
                thermal_diffusivity * (radial_conduction + cylindrical_term)
                + source_term
            )

        # Surface convection boundary condition
        T[n, -1] = (
            (conductivity / dr) * T[n, -2]
            + convection_coefficient_value * ambient_temp
        ) / (conductivity / dr + convection_coefficient_value)

    return T, r, t


# Crank-Nicolson solver
def crank_nicolson_solver(
    n_nodes,
    dt,
    t_end,
    radius=wire_radius,
    density=rho,
    heat_capacity=specific_heat,
    conductivity=thermal_conductivity,
    volumetric_heat_generation=heat_generation,
    initial_temp=initial_temperature,
    ambient_temp=ambient_temperature,
    convection_coefficient_value=convection_coefficient,
):
    thermal_diffusivity = conductivity / (density * heat_capacity)
    source_term = volumetric_heat_generation / (density * heat_capacity)

    r, t, dr = create_grid(n_nodes, dt, t_end, radius)
    n_time_steps = len(t)

    T = np.zeros((n_time_steps, n_nodes))
    T[0, :] = initial_temp

    # Coefficient matrix
    A = np.zeros((n_nodes, n_nodes))

    # Centreline symmetry boundary condition
    centre_coefficient = thermal_diffusivity * dt / dr**2

    A[0, 0] = 1 + 2 * centre_coefficient
    A[0, 1] = -2 * centre_coefficient

    # Interior nodes
    for j in range(1, n_nodes - 1):
        r_j = r[j]

        lower = (
            -thermal_diffusivity * dt / (2 * dr**2)
            + thermal_diffusivity * dt / (4 * r_j * dr)
        )
        diagonal = 1 + thermal_diffusivity * dt / dr**2
        upper = (
            -thermal_diffusivity * dt / (2 * dr**2)
            - thermal_diffusivity * dt / (4 * r_j * dr)
        )

        A[j, j - 1] = lower
        A[j, j] = diagonal
        A[j, j + 1] = upper

    # Surface convection boundary condition
    A[-1, -2] = conductivity / dr
    A[-1, -1] = -conductivity / dr - convection_coefficient_value

    for n in range(1, n_time_steps):
        previous_T = T[n - 1]
        b = np.zeros(n_nodes)

        # Centreline symmetry boundary condition
        b[0] = (
            previous_T[0]
            + 2
            * centre_coefficient
            * (previous_T[1] - previous_T[0])
            + source_term * dt
        )

        # Interior nodes
        for j in range(1, n_nodes - 1):
            radial_conduction = (
                previous_T[j + 1] - 2 * previous_T[j] + previous_T[j - 1]
            ) / dr**2

            cylindrical_term = (previous_T[j + 1] - previous_T[j - 1]) / (
                2 * r[j] * dr
            )

            b[j] = (
                previous_T[j]
                + thermal_diffusivity
                * dt
                / 2
                * (radial_conduction + cylindrical_term)
                + source_term * dt
            )

        # Surface convection boundary condition
        b[-1] = -convection_coefficient_value * ambient_temp

        T[n, :] = np.linalg.solve(A, b)

    return T, r, t


# Calculate Fourier number
def fourier_number(
    n_nodes,
    dt,
    radius=wire_radius,
    density=rho,
    heat_capacity=specific_heat,
    conductivity=thermal_conductivity,
):
    thermal_diffusivity = conductivity / (density * heat_capacity)
    dr = radius / (n_nodes - 1)

    return thermal_diffusivity * dt / dr**2