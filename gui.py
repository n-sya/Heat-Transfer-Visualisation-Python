import tkinter as tk
from tkinter import messagebox, ttk

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from calculations import crank_nicolson_solver, explicit_solver, fourier_number
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


class HeatTransferGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Transient Heat Conduction")
        self.root.geometry("1400x850")

        self.animation = None

        self.create_layout()

    # Create main GUI layout
    def create_layout(self):
        self.input_frame = ttk.Frame(self.root, padding=15)
        self.input_frame.pack(side=tk.LEFT, fill=tk.Y)

        self.results_frame = ttk.Frame(self.root, padding=10)
        self.results_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.create_input_panel()
        self.create_results_panel()

    # Create parameter input panel
    def create_input_panel(self):
        ttk.Label(
            self.input_frame,
            text="Simulation Inputs",
            font=("Segoe UI", 14, "bold"),
        ).grid(row=0, column=0, columnspan=2, pady=(0, 15))

        row = 1

        ttk.Label(
            self.input_frame,
            text="Geometry",
            font=("Segoe UI", 10, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(5, 5))

        row += 1
        self.radius_entry = self.create_entry(
            "Wire radius (m)",
            wire_radius,
            row,
        )

        row += 1

        ttk.Label(
            self.input_frame,
            text="Material Properties",
            font=("Segoe UI", 10, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(15, 5))

        row += 1
        self.conductivity_entry = self.create_entry(
            "Thermal conductivity (W/m K)",
            thermal_conductivity,
            row,
        )

        row += 1
        self.density_entry = self.create_entry(
            "Density (kg/m³)",
            rho,
            row,
        )

        row += 1
        self.heat_capacity_entry = self.create_entry(
            "Specific heat (J/kg K)",
            specific_heat,
            row,
        )

        row += 1

        ttk.Label(
            self.input_frame,
            text="Thermal Conditions",
            font=("Segoe UI", 10, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(15, 5))

        row += 1
        self.heat_generation_entry = self.create_entry(
            "Heat generation (W/m³)",
            heat_generation,
            row,
        )

        row += 1
        self.initial_temperature_entry = self.create_entry(
            "Initial temperature (K)",
            initial_temperature,
            row,
        )

        row += 1
        self.ambient_temperature_entry = self.create_entry(
            "Ambient temperature (K)",
            ambient_temperature,
            row,
        )

        row += 1
        self.convection_entry = self.create_entry(
            "Convection coefficient (W/m² K)",
            convection_coefficient,
            row,
        )

        row += 1

        ttk.Label(
            self.input_frame,
            text="Numerical Parameters",
            font=("Segoe UI", 10, "bold"),
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(15, 5))

        row += 1
        self.nodes_entry = self.create_entry(
            "Radial nodes",
            11,
            row,
        )

        row += 1
        self.dt_entry = self.create_entry(
            "Time step (s)",
            0.001,
            row,
        )

        row += 1
        self.duration_entry = self.create_entry(
            "Duration (s)",
            130,
            row,
        )

        row += 1

        ttk.Button(
            self.input_frame,
            text="Run Simulation",
            command=self.run_simulation,
        ).grid(
            row=row,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(20, 5),
        )

        row += 1

        self.status_label = ttk.Label(
            self.input_frame,
            text="Ready",
            wraplength=280,
        )
        self.status_label.grid(
            row=row,
            column=0,
            columnspan=2,
            sticky="w",
            pady=(10, 0),
        )

    # Create labelled numerical entry
    def create_entry(self, label, default_value, row):
        ttk.Label(
            self.input_frame,
            text=label,
        ).grid(
            row=row,
            column=0,
            sticky="w",
            padx=(0, 10),
            pady=3,
        )

        entry = ttk.Entry(
            self.input_frame,
            width=14,
        )
        entry.insert(0, str(default_value))
        entry.grid(
            row=row,
            column=1,
            sticky="ew",
            pady=3,
        )

        return entry

    # Create result tabs
    def create_results_panel(self):
        self.notebook = ttk.Notebook(self.results_frame)
        self.notebook.pack(
            fill=tk.BOTH,
            expand=True,
        )

        self.comparison_tab = ttk.Frame(self.notebook)
        self.distribution_tab = ttk.Frame(self.notebook)
        self.heatmap_tab = ttk.Frame(self.notebook)
        self.animation_tab = ttk.Frame(self.notebook)

        self.notebook.add(
            self.comparison_tab,
            text="Solver Comparison",
        )

        self.notebook.add(
            self.distribution_tab,
            text="Temperature Distribution",
        )

        self.notebook.add(
            self.heatmap_tab,
            text="Heatmap",
        )

        self.notebook.add(
            self.animation_tab,
            text="Animation",
        )

        self.create_placeholder(self.comparison_tab)
        self.create_placeholder(self.distribution_tab)
        self.create_placeholder(self.heatmap_tab)
        self.create_placeholder(self.animation_tab)

    # Create initial tab message
    def create_placeholder(self, parent):
        ttk.Label(
            parent,
            text="Run the simulation to display results.",
            font=("Segoe UI", 11),
        ).pack(expand=True)

    # Read and validate GUI inputs
    def get_inputs(self):
        radius = float(self.radius_entry.get())
        density = float(self.density_entry.get())
        heat_capacity = float(self.heat_capacity_entry.get())
        conductivity = float(self.conductivity_entry.get())
        volumetric_heat_generation = float(self.heat_generation_entry.get())
        initial_temp = float(self.initial_temperature_entry.get())
        ambient_temp = float(self.ambient_temperature_entry.get())
        convection_value = float(self.convection_entry.get())
        n_nodes = int(self.nodes_entry.get())
        dt = float(self.dt_entry.get())
        t_end = float(self.duration_entry.get())

        if radius <= 0:
            raise ValueError("Wire radius must be greater than zero.")

        if density <= 0:
            raise ValueError("Density must be greater than zero.")

        if heat_capacity <= 0:
            raise ValueError("Specific heat must be greater than zero.")

        if conductivity <= 0:
            raise ValueError("Thermal conductivity must be greater than zero.")

        if volumetric_heat_generation < 0:
            raise ValueError("Heat generation cannot be negative.")

        if convection_value < 0:
            raise ValueError("Convection coefficient cannot be negative.")

        if n_nodes < 3:
            raise ValueError("At least three radial nodes are required.")

        if dt <= 0:
            raise ValueError("Time step must be greater than zero.")

        if t_end <= 0:
            raise ValueError("Simulation duration must be greater than zero.")

        return {
            "radius": radius,
            "density": density,
            "heat_capacity": heat_capacity,
            "conductivity": conductivity,
            "volumetric_heat_generation": volumetric_heat_generation,
            "initial_temp": initial_temp,
            "ambient_temp": ambient_temp,
            "convection_coefficient_value": convection_value,
            "n_nodes": n_nodes,
            "dt": dt,
            "t_end": t_end,
        }

    # Run both numerical solvers
    def run_simulation(self):
        try:
            parameters = self.get_inputs()

        except ValueError as error:
            messagebox.showerror(
                "Invalid Input",
                str(error),
            )
            return

        self.status_label.config(text="Running simulation...")
        self.root.update_idletasks()

        solver_parameters = {
            key: value
            for key, value in parameters.items()
            if key not in {"n_nodes", "dt", "t_end"}
        }

        try:
            self.T_explicit, self.r, self.t = explicit_solver(
                parameters["n_nodes"],
                parameters["dt"],
                parameters["t_end"],
                **solver_parameters,
            )

            self.T_crank_nicolson, _, _ = crank_nicolson_solver(
                parameters["n_nodes"],
                parameters["dt"],
                parameters["t_end"],
                **solver_parameters,
            )

        except Exception as error:
            messagebox.showerror(
                "Simulation Error",
                str(error),
            )
            self.status_label.config(text="Simulation failed.")
            return

        self.parameters = parameters

        Fo = fourier_number(
            parameters["n_nodes"],
            parameters["dt"],
            radius=parameters["radius"],
            density=parameters["density"],
            heat_capacity=parameters["heat_capacity"],
            conductivity=parameters["conductivity"],
        )

        maximum_difference = np.max(
            np.abs(self.T_explicit[-1] - self.T_crank_nicolson[-1])
        )

        stability_status = "Stable" if Fo <= 0.5 else "Stability limit exceeded"

        self.status_label.config(
            text=(
                f"Simulation complete\n\n"
                f"Fourier number: {Fo:.4f}\n"
                f"Explicit scheme: {stability_status}\n\n"
                f"Explicit centre: {self.T_explicit[-1, 0]:.3f} K\n"
                f"Crank-Nicolson centre: "
                f"{self.T_crank_nicolson[-1, 0]:.3f} K\n\n"
                f"Maximum solver difference: "
                f"{maximum_difference:.6f} K"
            )
        )

        self.plot_solver_comparison()
        self.plot_temperature_distribution()
        self.plot_heatmaps()
        self.create_animation()

        # Open animation tab after simulation
        self.notebook.select(self.animation_tab)

    # Remove existing widgets from a results tab
    def clear_tab(self, tab):
        for widget in tab.winfo_children():
            widget.destroy()

    # Plot absolute difference between numerical methods
    def plot_solver_comparison(self):
        self.clear_tab(self.comparison_tab)

        temperature_difference = np.abs(
            self.T_explicit[-1] - self.T_crank_nicolson[-1]
        )

        figure, ax = plt.subplots(
            figsize=(8, 6),
            constrained_layout=True,
        )

        ax.plot(
            self.r * 1000,
            temperature_difference,
        )

        ax.set_xlabel("Radius (mm)")
        ax.set_ylabel("Absolute Difference (K)")
        ax.set_title("Difference Between Numerical Methods")
        ax.grid()

        canvas = FigureCanvasTkAgg(
            figure,
            master=self.comparison_tab,
        )
        canvas.draw()
        canvas.get_tk_widget().pack(
            fill=tk.BOTH,
            expand=True,
        )

    # Plot full temperature distribution from -R to R
    def plot_temperature_distribution(self):
        self.clear_tab(self.distribution_tab)

        explicit_full = np.concatenate(
            (
                self.T_explicit[-1][::-1],
                self.T_explicit[-1][1:],
            )
        )

        crank_nicolson_full = np.concatenate(
            (
                self.T_crank_nicolson[-1][::-1],
                self.T_crank_nicolson[-1][1:],
            )
        )

        r_full = np.concatenate(
            (
                -self.r[::-1],
                self.r[1:],
            )
        )

        figure, ax = plt.subplots(
            figsize=(8, 6),
            constrained_layout=True,
        )

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

        ax.set_xlabel("Radius (mm)")
        ax.set_ylabel("Temperature (K)")
        ax.set_title("Final Temperature Distribution Across Diameter")
        ax.grid()
        ax.legend()

        canvas = FigureCanvasTkAgg(
            figure,
            master=self.distribution_tab,
        )
        canvas.draw()
        canvas.get_tk_widget().pack(
            fill=tk.BOTH,
            expand=True,
        )

    # Create 2D radial temperature field
    def create_temperature_field(self, temperature_profile):
        radius = self.parameters["radius"]
        points = 200

        x = np.linspace(-radius, radius, points)
        y = np.linspace(-radius, radius, points)

        X, Y = np.meshgrid(x, y)

        radial_distance = np.sqrt(X**2 + Y**2)

        temperature_field = np.interp(
            radial_distance,
            self.r,
            temperature_profile,
        )

        temperature_field[radial_distance > radius] = np.nan

        return temperature_field

    # Create common temperature limits
    def get_temperature_limits(self):
        minimum_temperature = min(
            np.nanmin(self.T_explicit),
            np.nanmin(self.T_crank_nicolson),
        )

        maximum_temperature = max(
            np.nanmax(self.T_explicit),
            np.nanmax(self.T_crank_nicolson),
        )

        if np.isclose(minimum_temperature, maximum_temperature):
            maximum_temperature = minimum_temperature + 1e-6

        return minimum_temperature, maximum_temperature

    # Plot final Explicit and Crank-Nicolson heatmaps
    def plot_heatmaps(self):
        self.clear_tab(self.heatmap_tab)

        explicit_field = self.create_temperature_field(
            self.T_explicit[-1]
        )

        crank_nicolson_field = self.create_temperature_field(
            self.T_crank_nicolson[-1]
        )

        radius_mm = self.parameters["radius"] * 1000

        minimum_temperature = min(
            np.nanmin(explicit_field),
            np.nanmin(crank_nicolson_field),
        )

        maximum_temperature = max(
            np.nanmax(explicit_field),
            np.nanmax(crank_nicolson_field),
        )

        if np.isclose(minimum_temperature, maximum_temperature):
            maximum_temperature = minimum_temperature + 1e-6

        levels = np.linspace(
            minimum_temperature,
            maximum_temperature,
            12,
        )

        points = explicit_field.shape[0]

        x = np.linspace(
            -radius_mm,
            radius_mm,
            points,
        )

        y = np.linspace(
            -radius_mm,
            radius_mm,
            points,
        )

        X, Y = np.meshgrid(x, y)

        figure = plt.figure(
            figsize=(11, 5),
            constrained_layout=True,
        )

        grid = figure.add_gridspec(
            1,
            3,
            width_ratios=[1, 1, 0.05],
        )

        explicit_ax = figure.add_subplot(grid[0, 0])
        crank_ax = figure.add_subplot(grid[0, 1])
        colorbar_ax = figure.add_subplot(grid[0, 2])

        explicit_contour = explicit_ax.contourf(
            X,
            Y,
            explicit_field,
            levels=levels,
        )

        crank_ax.contourf(
            X,
            Y,
            crank_nicolson_field,
            levels=levels,
        )

        explicit_ax.set_title("Explicit")
        crank_ax.set_title("Crank-Nicolson")

        for ax in (explicit_ax, crank_ax):
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            ax.set_aspect("equal")

        colorbar = figure.colorbar(
            explicit_contour,
            cax=colorbar_ax,
        )

        colorbar.set_label("Temperature (K)")
        colorbar.formatter = ticker.FormatStrFormatter("%.3f")
        colorbar.update_ticks()

        canvas = FigureCanvasTkAgg(
            figure,
            master=self.heatmap_tab,
        )
        canvas.draw()
        canvas.get_tk_widget().pack(
            fill=tk.BOTH,
            expand=True,
        )

    # Create animated Explicit and Crank-Nicolson heatmaps
    def create_animation(self):
        self.clear_tab(self.animation_tab)

        radius_mm = self.parameters["radius"] * 1000

        frame_indices = np.linspace(
            0,
            len(self.t) - 1,
            min(150, len(self.t)),
            dtype=int,
        )

        minimum_temperature, maximum_temperature = (
            self.get_temperature_limits()
        )

        initial_explicit_field = self.create_temperature_field(
            self.T_explicit[frame_indices[0]]
        )

        initial_crank_field = self.create_temperature_field(
            self.T_crank_nicolson[frame_indices[0]]
        )

        figure = plt.figure(
            figsize=(11, 5),
            constrained_layout=True,
        )

        grid = figure.add_gridspec(
            1,
            3,
            width_ratios=[1, 1, 0.05],
        )

        explicit_ax = figure.add_subplot(grid[0, 0])
        crank_ax = figure.add_subplot(grid[0, 1])
        colorbar_ax = figure.add_subplot(grid[0, 2])

        explicit_image = explicit_ax.imshow(
            initial_explicit_field,
            origin="lower",
            extent=[
                -radius_mm,
                radius_mm,
                -radius_mm,
                radius_mm,
            ],
            vmin=minimum_temperature,
            vmax=maximum_temperature,
        )

        crank_image = crank_ax.imshow(
            initial_crank_field,
            origin="lower",
            extent=[
                -radius_mm,
                radius_mm,
                -radius_mm,
                radius_mm,
            ],
            vmin=minimum_temperature,
            vmax=maximum_temperature,
        )

        explicit_ax.set_title("Explicit")
        crank_ax.set_title("Crank-Nicolson")

        for ax in (explicit_ax, crank_ax):
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            ax.set_aspect("equal")

        colorbar = figure.colorbar(
            explicit_image,
            cax=colorbar_ax,
        )

        colorbar.set_label("Temperature (K)")
        colorbar.formatter = ticker.FormatStrFormatter("%.2f")
        colorbar.update_ticks()

        time_text = figure.suptitle(
            f"Time = {self.t[frame_indices[0]]:.2f} s"
        )

        def update(frame_number):
            index = frame_indices[frame_number]

            explicit_field = self.create_temperature_field(
                self.T_explicit[index]
            )

            crank_field = self.create_temperature_field(
                self.T_crank_nicolson[index]
            )

            explicit_image.set_data(explicit_field)
            crank_image.set_data(crank_field)

            time_text.set_text(
                f"Time = {self.t[index]:.2f} s"
            )

            return (
                explicit_image,
                crank_image,
                time_text,
            )

        self.animation = animation.FuncAnimation(
            figure,
            update,
            frames=len(frame_indices),
            interval=60,
            repeat=True,
            blit=False,
        )

        canvas = FigureCanvasTkAgg(
            figure,
            master=self.animation_tab,
        )
        canvas.draw()
        canvas.get_tk_widget().pack(
            fill=tk.BOTH,
            expand=True,
        )

        controls = ttk.Frame(self.animation_tab)
        controls.pack(pady=5)

        ttk.Button(
            controls,
            text="Play",
            command=self.play_animation,
        ).pack(
            side=tk.LEFT,
            padx=5,
        )

        ttk.Button(
            controls,
            text="Pause",
            command=self.pause_animation,
        ).pack(
            side=tk.LEFT,
            padx=5,
        )

        ttk.Button(
            controls,
            text="Restart",
            command=self.restart_animation,
        ).pack(
            side=tk.LEFT,
            padx=5,
        )

    # Play animation
    def play_animation(self):
        if self.animation is not None:
            self.animation.event_source.start()

    # Pause animation
    def pause_animation(self):
        if self.animation is not None:
            self.animation.event_source.stop()

    # Restart animation
    def restart_animation(self):
        if self.animation is not None:
            self.animation.frame_seq = self.animation.new_frame_seq()
            self.animation.event_source.start()


# Launch application
if __name__ == "__main__":
    root = tk.Tk()
    app = HeatTransferGUI(root)
    root.mainloop()