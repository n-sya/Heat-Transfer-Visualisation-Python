import tkinter as tk
import unittest

import numpy as np

from calculations import (
    crank_nicolson_solver,
    create_grid,
    explicit_solver,
    fourier_number,
)
from gui import HeatTransferGUI
from variables import (
    ambient_temperature,
    initial_temperature,
    rho,
    specific_heat,
    thermal_conductivity,
    wire_radius,
)


class TestCalculations(unittest.TestCase):
    def test_create_grid_reaches_wire_surface(self):
        n_nodes = 11
        dt = 0.001
        t_end = 1.0

        r, t, dr = create_grid(n_nodes, dt, t_end)

        self.assertEqual(len(r), n_nodes)
        self.assertAlmostEqual(r[0], 0.0)
        self.assertAlmostEqual(r[-1], wire_radius)
        self.assertAlmostEqual(dr, wire_radius / (n_nodes - 1))
        self.assertAlmostEqual(t[0], 0.0)
        self.assertAlmostEqual(t[-1], t_end)

    def test_radial_grid_is_uniform(self):
        r, _, dr = create_grid(11, 0.001, 1.0)

        radial_spacing = np.diff(r)

        self.assertTrue(np.allclose(radial_spacing, dr))

    def test_fourier_number_matches_manual_calculation(self):
        n_nodes = 11
        dt = 0.001

        thermal_diffusivity = thermal_conductivity / (rho * specific_heat)
        dr = wire_radius / (n_nodes - 1)

        expected_fourier_number = thermal_diffusivity * dt / dr**2
        calculated_fourier_number = fourier_number(n_nodes, dt)

        self.assertAlmostEqual(
            calculated_fourier_number,
            expected_fourier_number,
        )

    def test_default_fourier_number_is_stable(self):
        Fo = fourier_number(11, 0.001)

        self.assertLessEqual(Fo, 0.5)

    def test_larger_time_step_increases_fourier_number(self):
        smaller_time_step = fourier_number(11, 0.001)
        larger_time_step = fourier_number(11, 0.002)

        self.assertGreater(larger_time_step, smaller_time_step)

    def test_explicit_solver_output_shape(self):
        T, r, t = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
        )

        self.assertEqual(T.shape, (len(t), 11))
        self.assertEqual(len(r), 11)

    def test_explicit_initial_condition(self):
        T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=0.01,
        )

        self.assertTrue(
            np.allclose(
                T[0],
                initial_temperature,
            )
        )

    def test_explicit_temperatures_are_finite(self):
        T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
        )

        self.assertTrue(np.all(np.isfinite(T)))

    def test_explicit_temperature_rises(self):
        T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
        )

        self.assertGreater(
            T[-1, 0],
            initial_temperature,
        )

    def test_explicit_centre_hotter_than_surface(self):
        T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=5.0,
        )

        self.assertGreaterEqual(
            T[-1, 0],
            T[-1, -1],
        )

    def test_crank_nicolson_solver_output_shape(self):
        T, r, t = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=1.0,
        )

        self.assertEqual(T.shape, (len(t), 11))
        self.assertEqual(len(r), 11)

    def test_crank_nicolson_initial_condition(self):
        T, _, _ = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=0.1,
        )

        self.assertTrue(
            np.allclose(
                T[0],
                initial_temperature,
            )
        )

    def test_crank_nicolson_temperatures_are_finite(self):
        T, _, _ = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=1.0,
        )

        self.assertTrue(np.all(np.isfinite(T)))

    def test_crank_nicolson_centreline_symmetry(self):
        T, _, _ = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=1.0,
        )

        self.assertTrue(
            np.allclose(
                T[:, 0],
                T[:, 1],
            )
        )

    def test_crank_nicolson_temperature_rises(self):
        T, _, _ = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=1.0,
        )

        self.assertGreater(
            T[-1, 0],
            initial_temperature,
        )

    def test_crank_nicolson_centre_hotter_than_surface(self):
        T, _, _ = crank_nicolson_solver(
            n_nodes=11,
            dt=0.01,
            t_end=5.0,
        )

        self.assertGreaterEqual(
            T[-1, 0],
            T[-1, -1],
        )


class TestFunctionalBehaviour(unittest.TestCase):
    def test_explicit_and_crank_nicolson_are_similar(self):
        n_nodes = 11
        dt = 0.001
        t_end = 5.0

        T_explicit, r_explicit, t_explicit = explicit_solver(
            n_nodes,
            dt,
            t_end,
        )

        T_crank_nicolson, r_crank, t_crank = crank_nicolson_solver(
            n_nodes,
            dt,
            t_end,
        )

        self.assertTrue(
            np.allclose(
                r_explicit,
                r_crank,
            )
        )

        self.assertTrue(
            np.allclose(
                t_explicit,
                t_crank,
            )
        )

        maximum_difference = np.max(
            np.abs(
                T_explicit[-1]
                - T_crank_nicolson[-1]
            )
        )

        self.assertLess(
            maximum_difference,
            0.1,
        )

    def test_temperature_distribution_is_symmetric(self):
        T, r, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
        )

        final_temperature = T[-1]

        full_temperature = np.concatenate(
            (
                final_temperature[::-1],
                final_temperature[1:],
            )
        )

        full_radius = np.concatenate(
            (
                -r[::-1],
                r[1:],
            )
        )

        midpoint = len(full_radius) // 2

        left_temperature = full_temperature[:midpoint]
        right_temperature = full_temperature[midpoint + 1 :][::-1]

        self.assertTrue(
            np.allclose(
                left_temperature,
                right_temperature,
            )
        )

    def test_modified_conductivity_changes_solution(self):
        default_T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
        )

        modified_T, _, _ = explicit_solver(
            n_nodes=11,
            dt=0.001,
            t_end=1.0,
            conductivity=thermal_conductivity * 2,
        )

        self.assertFalse(
            np.allclose(
                default_T[-1],
                modified_T[-1],
            )
        )

    def test_complete_simulation_workflow(self):
        n_nodes = 11
        dt = 0.001
        t_end = 1.0

        T_explicit, r_explicit, t_explicit = explicit_solver(
            n_nodes,
            dt,
            t_end,
        )

        T_crank_nicolson, r_crank, t_crank = crank_nicolson_solver(
            n_nodes,
            dt,
            t_end,
        )

        Fo = fourier_number(
            n_nodes,
            dt,
        )

        difference = np.abs(
            T_explicit[-1]
            - T_crank_nicolson[-1]
        )

        self.assertEqual(
            T_explicit.shape,
            T_crank_nicolson.shape,
        )

        self.assertTrue(
            np.allclose(
                r_explicit,
                r_crank,
            )
        )

        self.assertTrue(
            np.allclose(
                t_explicit,
                t_crank,
            )
        )

        self.assertLessEqual(Fo, 0.5)
        self.assertTrue(np.all(np.isfinite(T_explicit)))
        self.assertTrue(np.all(np.isfinite(T_crank_nicolson)))
        self.assertLess(np.max(difference), 0.1)


class TestGUI(unittest.TestCase):
    def setUp(self):
        self.root = tk.Tk()
        self.root.withdraw()
        self.gui = HeatTransferGUI(self.root)

    def tearDown(self):
        self.root.destroy()

    def set_entry(self, entry, value):
        entry.delete(0, tk.END)
        entry.insert(0, str(value))

    def test_default_gui_inputs(self):
        parameters = self.gui.get_inputs()

        self.assertAlmostEqual(
            parameters["radius"],
            wire_radius,
        )
        self.assertAlmostEqual(
            parameters["density"],
            rho,
        )
        self.assertAlmostEqual(
            parameters["heat_capacity"],
            specific_heat,
        )
        self.assertAlmostEqual(
            parameters["conductivity"],
            thermal_conductivity,
        )
        self.assertAlmostEqual(
            parameters["initial_temp"],
            initial_temperature,
        )
        self.assertAlmostEqual(
            parameters["ambient_temp"],
            ambient_temperature,
        )

        self.assertEqual(
            parameters["n_nodes"],
            11,
        )
        self.assertAlmostEqual(
            parameters["dt"],
            0.001,
        )
        self.assertAlmostEqual(
            parameters["t_end"],
            130.0,
        )

    def test_gui_rejects_negative_radius(self):
        self.set_entry(
            self.gui.radius_entry,
            -0.001,
        )

        with self.assertRaises(ValueError):
            self.gui.get_inputs()

    def test_gui_rejects_invalid_number_of_nodes(self):
        self.set_entry(
            self.gui.nodes_entry,
            2,
        )

        with self.assertRaises(ValueError):
            self.gui.get_inputs()

    def test_gui_rejects_zero_time_step(self):
        self.set_entry(
            self.gui.dt_entry,
            0,
        )

        with self.assertRaises(ValueError):
            self.gui.get_inputs()

    def test_gui_rejects_negative_heat_generation(self):
        self.set_entry(
            self.gui.heat_generation_entry,
            -1,
        )

        with self.assertRaises(ValueError):
            self.gui.get_inputs()

    def test_temperature_field_shape(self):
        self.gui.parameters = {
            "radius": wire_radius,
        }

        self.gui.r = np.linspace(
            0,
            wire_radius,
            11,
        )

        temperature_profile = np.linspace(
            310,
            300,
            11,
        )

        temperature_field = self.gui.create_temperature_field(
            temperature_profile
        )

        self.assertEqual(
            temperature_field.shape,
            (200, 200),
        )

    def test_temperature_field_masks_outside_wire(self):
        self.gui.parameters = {
            "radius": wire_radius,
        }

        self.gui.r = np.linspace(
            0,
            wire_radius,
            11,
        )

        temperature_profile = np.linspace(
            310,
            300,
            11,
        )

        temperature_field = self.gui.create_temperature_field(
            temperature_profile
        )

        self.assertTrue(
            np.isnan(
                temperature_field[0, 0]
            )
        )
        self.assertTrue(
            np.isnan(
                temperature_field[0, -1]
            )
        )
        self.assertTrue(
            np.isnan(
                temperature_field[-1, 0]
            )
        )
        self.assertTrue(
            np.isnan(
                temperature_field[-1, -1]
            )
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)