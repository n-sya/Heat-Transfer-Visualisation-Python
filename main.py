from calculations import (
    crank_nicolson_solver,
    explicit_solver,
    fourier_number,
)
from visualisation import (
    animate_temperature,
    plot_solver_comparison,
    plot_temperature_heatmap,
    plot_temperature_profile,
)

# Numerical parameters
n_nodes = 11
dt = 0.001
t_end = 130

# Run solvers
T_explicit, r, t = explicit_solver(n_nodes, dt, t_end)
T_crank_nicolson, _, _ = crank_nicolson_solver(n_nodes, dt, t_end)

# Numerical stability
Fo = fourier_number(n_nodes, dt)

print(f"Fourier number: {Fo:.4f}")

# Generate explicit solver outputs
plot_temperature_profile(T_explicit, r, "Explicit")
plot_temperature_heatmap(T_explicit, r, "Explicit")
animate_temperature(T_explicit, r, t, "Explicit")

# Generate Crank-Nicolson solver outputs
plot_temperature_profile(T_crank_nicolson, r, "Crank-Nicolson")
plot_temperature_heatmap(T_crank_nicolson, r, "Crank-Nicolson")
animate_temperature(T_crank_nicolson, r, t, "Crank-Nicolson")

# Compare solvers
plot_solver_comparison(T_explicit, T_crank_nicolson, r)