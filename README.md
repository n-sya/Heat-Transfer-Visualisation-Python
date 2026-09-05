# Transient Heat Conduction in a Cylindrical Resistor

A Python simulation of transient radial heat conduction in a cylindrical resistor with internal heat generation and convective cooling. The model implements and compares Explicit Finite Difference and Crank-Nicolson numerical methods, with an interactive GUI for modifying simulation parameters and visualising the resulting temperature field.

## Overview

This project models the transient temperature distribution within a cylindrical resistor subjected to uniform volumetric heat generation. Heat is conducted radially through the resistor and rejected to the surroundings through convection at its outer surface.

The transient heat equation is solved numerically using two methods:

- Explicit Finite Difference Method
- Crank-Nicolson Method

The implementation allows the numerical behaviour of the two methods to be compared while demonstrating the stability limitations associated with the Explicit scheme.

An interactive Tkinter GUI provides control over the physical, thermal, and numerical parameters of the simulation and displays the resulting temperature distributions and heatmaps.

## Features

- Transient one-dimensional radial heat conduction
- Cylindrical coordinate formulation
- Uniform volumetric heat generation
- Centreline symmetry boundary condition
- Convective outer boundary condition
- Explicit Finite Difference solver
- Crank-Nicolson solver
- Explicit stability checking using the Fourier number
- Comparison between numerical methods
- Interactive simulation parameters
- Temperature distribution across the resistor diameter
- Two-dimensional temperature heatmaps
- Animated temperature evolution
- Automated unit and functional testing
- JSON-driven simulation test cases

## Mathematical Model

For axisymmetric radial conduction with constant material properties, the governing equation is

$$
\rho c_p \frac{\partial T}{\partial t} = k \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) + \dot{q}
$$

where:

- $T$ is temperature
- $t$ is time
- $r$ is radial position
- $\rho$ is density
- $c_p$ is specific heat capacity
- $k$ is thermal conductivity
- $\dot{q}$ is volumetric heat generation

The thermal diffusivity is defined as

$$
\alpha = \frac{k}{\rho c_p}
$$

giving

$$
\frac{\partial T}{\partial t}=\alpha\left(\frac{\partial^2 T}{\partial r^2}+\frac{1}{r}\frac{\partial T}{\partial r}\right)+\frac{\dot{q}}{\rho c_p}
$$

### Initial Condition

The resistor is initially assumed to have a uniform temperature:

$$
T(r,0)=T_0
$$

### Centreline Boundary Condition

Symmetry at the centre of the cylinder requires

$$
\left.\frac{\partial T}{\partial r}\right|_{r=0}=0
$$

### Surface Boundary Condition

At the outer surface, heat is transferred to the surroundings through convection:

$$
-k\left.\frac{\partial T}{\partial r}\right|_{r=R}=h(T_s-T_\infty)
$$

where $h$ is the convection coefficient and $T_\infty$ is the ambient temperature.

A full derivation of the numerical schemes and boundary-condition discretisations is provided in [DERIVATIONS.md](DERIVATIONS.md).

## Numerical Methods

### Explicit Finite Difference Method

The Explicit method calculates the temperature at the next time step directly from temperatures at the current time step.

This makes the method straightforward to implement but conditionally stable. The stability of the implemented scheme is assessed using the Fourier number

$$
Fo = \frac{\alpha \Delta t}{\Delta r^2}
$$

with the stability requirement

$$
Fo \leq 0.25
$$

The centreline node imposes the most restrictive stability condition for the cylindrical formulation used in this model.

The GUI evaluates this condition before running the simulation. If the selected parameters exceed the stability limit, the simulation is stopped and the maximum allowable time step is reported.

### Crank-Nicolson Method

The Crank-Nicolson method evaluates the spatial conduction terms using the average of the current and next time levels.

This produces a system of simultaneous linear equations of the form

$$
A\mathbf{T}^{n+1}=\mathbf{b}
$$

which is solved using NumPy.

For the linear diffusion problem considered here, the Crank-Nicolson method is unconditionally stable, although the choice of time step continues to affect numerical accuracy.

## Interactive GUI

The graphical interface allows the following parameters to be modified:

**Geometry**

- Wire radius

**Material properties**

- Thermal conductivity
- Density
- Specific heat capacity

**Thermal conditions**

- Volumetric heat generation
- Initial temperature
- Ambient temperature
- Convection coefficient

**Numerical parameters**

- Number of radial nodes
- Time step
- Simulation duration

The results are displayed across four tabs:

1. **Solver Comparison**  
   Displays the absolute difference between the Explicit and Crank-Nicolson solutions.

2. **Temperature Distribution**  
   Displays the final temperature distribution across the complete resistor diameter.

3. **Heatmap**  
   Displays the final two-dimensional temperature fields predicted by both numerical methods.

4. **Animation**  
   Displays the transient evolution of both temperature fields with play, pause, and restart controls.

## Results

The simulation produces temperature profiles, heatmaps, animations, and a comparison between the two numerical methods.

Example outputs are stored in the `outputs/` directory.

The two numerical methods can be compared directly to assess agreement between the independently implemented discretisations.

## Repository Structure

```text
Heat-Transfer-Visualisation-Python/
├── archive/
│   ├── crank_nicolson_solver.py
│   └── explicit_solver.py
├── outputs/
├── calculations.py
├── variables.py
├── visualisation.py
├── gui.py
├── main.py
├── test.py
├── test_cases.json
├── DERIVATIONS.md
├── README.md
└── .gitignore
```

### `calculations.py`

Contains the numerical grid generation, Explicit solver, Crank-Nicolson solver, and Fourier number calculation.

### `variables.py`

Contains the default physical and thermal parameters used by the model.

### `visualisation.py`

Generates temperature profiles, heatmaps, animations, and solver-comparison plots.

### `gui.py`

Provides the interactive Tkinter interface for modifying simulation parameters and visualising results.

### `main.py`

Runs the simulation using the default parameters and generates the standard output files.

### `test.py`

Contains unit and functional tests for the numerical calculations and GUI behaviour.

### `test_cases.json`

Contains simulation configurations used for data-driven functional testing.

### `DERIVATIONS.md`

Contains the mathematical derivation of the finite difference equations implemented by the numerical solvers.

### `archive/`

Preserves the original solver implementations developed before the project was restructured.

## Running the Simulation

Clone the repository and move into the project directory.

Install the required Python packages:

```bash
pip install numpy matplotlib
```

### Run the GUI

```bash
python gui.py
```

### Generate Standard Outputs

```bash
python main.py
```

Generated plots and animations are stored in the `outputs/` directory.

## Testing

The project uses Python's built-in `unittest` framework.

Run the complete test suite using:

```bash
python test.py
```

The tests assess:

- Grid generation
- Fourier number calculation
- Explicit solver behaviour
- Crank-Nicolson solver behaviour
- Initial conditions
- Boundary behaviour
- Numerical solution dimensions
- Finite numerical outputs
- Agreement between the numerical methods
- GUI input handling
- Temperature-field generation
- JSON-defined simulation configurations

The JSON test cases allow different physical and numerical configurations to be evaluated without duplicating test code.

## Further Documentation

For the complete mathematical development of the numerical model, including the treatment of the centreline and convective boundary conditions, see:

**[Numerical Derivations](DERIVATIONS.md)**