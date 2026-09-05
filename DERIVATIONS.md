# Numerical Derivations

This document presents the mathematical formulation and finite difference discretisation used to model transient radial heat conduction in a cylindrical resistor.

The model considers one-dimensional radial heat conduction with uniform volumetric heat generation and convective heat transfer from the outer surface.

## 1. Governing Equation

For transient heat conduction in cylindrical coordinates, assuming axisymmetry and no axial or angular temperature variation, the governing equation is

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

Defining thermal diffusivity as

$$
\alpha = \frac{k}{\rho c_p}
$$

and the heat-generation source term as

$$
S = \frac{\dot{q}}{\rho c_p}
$$

the governing equation becomes

$$
\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) + S
$$

## 2. Initial Condition

The resistor is initially assumed to have a uniform temperature $T_0$:

$$
T(r,0) = T_0
$$

Therefore, at the initial time level,

$$
T_j^0 = T_0
$$

for every radial node $j$.

## 3. Boundary Conditions

Two boundary conditions are required because the spatial governing equation is second order.

### 3.1 Centreline Symmetry

At the centre of the cylinder, radial symmetry requires the temperature gradient to vanish:

$$
\left. \frac{\partial T}{\partial r} \right|_{r=0} = 0
$$

No heat therefore crosses the centreline.

### 3.2 Surface Convection

At the outer radius $R$, conductive heat transfer through the resistor surface is balanced by convection to the surroundings:

$$
-k \left. \frac{\partial T}{\partial r} \right|_{r=R} = h(T_s-T_\infty)
$$

where:

- $h$ is the convection coefficient
- $T_s$ is the resistor surface temperature
- $T_\infty$ is the ambient temperature

## 4. Numerical Grid

The radial domain

$$
0 \leq r \leq R
$$

is divided into $N$ equally spaced nodes.

The radial spacing is

$$
\Delta r = \frac{R}{N-1}
$$

and the position of node $j$ is

$$
r_j = j\Delta r
$$

where

$$
j = 0,1,\ldots,N-1
$$

The temporal domain is discretised using a constant time step $\Delta t$:

$$
t^n = n\Delta t
$$

The numerical temperature at radial node $j$ and time level $n$ is denoted

$$
T_j^n
$$

## 5. Spatial Finite Difference Approximations

For an interior node, the second radial derivative is approximated using a central difference:

$$
\left. \frac{\partial^2 T}{\partial r^2} \right|_j \approx \frac{T_{j+1}-2T_j+T_{j-1}}{\Delta r^2}
$$

The first radial derivative is also approximated using a central difference:

$$
\left. \frac{\partial T}{\partial r} \right|_j \approx \frac{T_{j+1}-T_{j-1}}{2\Delta r}
$$

Therefore, the cylindrical spatial operator becomes

$$
\left. \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) \right|_j \approx \frac{T_{j+1}-2T_j+T_{j-1}}{\Delta r^2} + \frac{1}{r_j}\frac{T_{j+1}-T_{j-1}}{2\Delta r}
$$

This expression is used for all interior nodes.

## 6. Explicit Finite Difference Method

The Explicit method approximates the temporal derivative using a forward difference:

$$
\frac{\partial T}{\partial t} \approx \frac{T_j^{n+1}-T_j^n}{\Delta t}
$$

Substitution into the governing equation gives

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t} = \alpha \left[ \frac{T_{j+1}^n-2T_j^n+T_{j-1}^n}{\Delta r^2} + \frac{T_{j+1}^n-T_{j-1}^n}{2r_j\Delta r} \right] + S
$$

Solving for the unknown temperature at the next time level gives

$$
T_j^{n+1} = T_j^n + \alpha\Delta t \left[ \frac{T_{j+1}^n-2T_j^n+T_{j-1}^n}{\Delta r^2} + \frac{T_{j+1}^n-T_{j-1}^n}{2r_j\Delta r} \right] + S\Delta t
$$

This is the Explicit update equation used for the interior nodes.

### 6.1 Fourier Number

The Fourier number is defined as

$$
Fo = \frac{\alpha\Delta t}{\Delta r^2}
$$

Using this definition, the Explicit interior equation can also be written as

$$
T_j^{n+1} = T_j^n + Fo(T_{j+1}^n-2T_j^n+T_{j-1}^n) + \frac{Fo\Delta r}{2r_j}(T_{j+1}^n-T_{j-1}^n) + S\Delta t
$$

Collecting the temperature coefficients gives

$$
T_j^{n+1} = Fo\left(1-\frac{\Delta r}{2r_j}\right)T_{j-1}^n + (1-2Fo)T_j^n + Fo\left(1+\frac{\Delta r}{2r_j}\right)T_{j+1}^n + S\Delta t
$$

## 7. Explicit Centreline Treatment

The cylindrical governing equation contains the term

$$
\frac{1}{r}\frac{\partial T}{\partial r}
$$

which appears singular at

$$
r=0
$$

and therefore requires separate treatment.

At the centreline, symmetry gives

$$
\left. \frac{\partial T}{\partial r} \right|_{r=0} = 0
$$

Consider the radial part of the cylindrical Laplacian:

$$
\frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r}
$$

Using the centreline limit,

$$
\lim_{r\rightarrow0}\frac{1}{r}\frac{\partial T}{\partial r} = \left. \frac{\partial^2 T}{\partial r^2} \right|_{r=0}
$$

so that

$$
\left. \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) \right|_{r=0} = 2\left. \frac{\partial^2 T}{\partial r^2} \right|_{r=0}
$$

Using a symmetric ghost node,

$$
T_{-1}^n = T_1^n
$$

the centreline curvature is approximated by

$$
\left. \frac{\partial^2 T}{\partial r^2} \right|_{r=0} \approx \frac{T_1^n-2T_0^n+T_{-1}^n}{\Delta r^2} = \frac{2(T_1^n-T_0^n)}{\Delta r^2}
$$

Therefore, the complete cylindrical spatial operator at the centreline becomes

$$
\left. \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) \right|_{r=0} \approx \frac{4(T_1^n-T_0^n)}{\Delta r^2}
$$

The Explicit centreline update is therefore

$$
T_0^{n+1} = T_0^n + \Delta t \left[ 4\alpha\frac{T_1^n-T_0^n}{\Delta r^2} + S \right]
$$

Using the Fourier number,

$$
T_0^{n+1} = T_0^n + 4Fo(T_1^n-T_0^n) + S\Delta t
$$

or

$$
T_0^{n+1} = (1-4Fo)T_0^n + 4FoT_1^n + S\Delta t
$$

This is the centreline equation implemented by the Explicit solver.

## 8. Explicit Surface Convection

At the outer surface,

$$
-k\left. \frac{\partial T}{\partial r} \right|_{r=R} = h(T_s-T_\infty)
$$

The radial temperature gradient is approximated using a backward difference:

$$
\left. \frac{\partial T}{\partial r} \right|_{r=R} \approx \frac{T_s-T_{N-2}}{\Delta r}
$$

Substitution gives

$$
-k\frac{T_s-T_{N-2}}{\Delta r} = h(T_s-T_\infty)
$$

Expanding,

$$
-\frac{k}{\Delta r}T_s + \frac{k}{\Delta r}T_{N-2} = hT_s-hT_\infty
$$

Rearranging,

$$
\left(\frac{k}{\Delta r}+h\right)T_s = \frac{k}{\Delta r}T_{N-2}+hT_\infty
$$

Therefore,

$$
T_s = \frac{\frac{k}{\Delta r}T_{N-2}+hT_\infty}{\frac{k}{\Delta r}+h}
$$

This expression is used to calculate the outer-surface temperature in the Explicit solver.

## 9. Explicit Stability

The Explicit method is conditionally stable.

For the interior radial nodes, the Explicit formulation gives the stability restriction

$$
Fo \leq 0.5
$$

However, the centreline update is

$$
T_0^{n+1} = (1-4Fo)T_0^n + 4FoT_1^n + S\Delta t
$$

For the coefficient of $T_0^n$ to remain non-negative,

$$
1-4Fo \geq 0
$$

and therefore

$$
Fo \leq 0.25
$$

The centreline condition is more restrictive than the interior-node condition and therefore determines the stability limit used by the Explicit solver.

Since

$$
Fo = \frac{\alpha\Delta t}{\Delta r^2}
$$

the stability requirement can be written as

$$
\frac{\alpha\Delta t}{\Delta r^2} \leq 0.25
$$

giving

$$
\Delta t \leq \frac{0.25\Delta r^2}{\alpha}
$$

Therefore,

$$
\boxed{\Delta t_{\max} = \frac{0.25\Delta r^2}{\alpha}}
$$

The GUI evaluates the Fourier number before running either numerical solver. If

$$
Fo > 0.25
$$

the simulation is stopped and the user is informed of the maximum permissible time step.

This prevents an unstable Explicit solution from being presented alongside the Crank-Nicolson solution.

## 10. Crank-Nicolson Method

The Crank-Nicolson method evaluates the spatial operator using the average between time levels $n$ and $n+1$.

Define the spatial operator

$$
L(T_j) = \frac{T_{j+1}-2T_j+T_{j-1}}{\Delta r^2} + \frac{T_{j+1}-T_{j-1}}{2r_j\Delta r}
$$

The governing equation can then be written as

$$
\frac{\partial T}{\partial t} = \alpha L(T) + S
$$

The Crank-Nicolson discretisation is

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t} = \frac{\alpha}{2}\left[L(T_j^{n+1})+L(T_j^n)\right]+S
$$

Multiplying by $\Delta t$ gives

$$
T_j^{n+1}-T_j^n = \frac{\alpha\Delta t}{2}L(T_j^{n+1}) + \frac{\alpha\Delta t}{2}L(T_j^n) + S\Delta t
$$

The unknown $n+1$ terms are moved to the left-hand side:

$$
T_j^{n+1} - \frac{\alpha\Delta t}{2}L(T_j^{n+1}) = T_j^n + \frac{\alpha\Delta t}{2}L(T_j^n) + S\Delta t
$$

## 11. Crank-Nicolson Interior Coefficients

Expanding the unknown spatial operator gives

$$
T_j^{n+1} - \frac{\alpha\Delta t}{2}\left[\frac{T_{j+1}^{n+1}-2T_j^{n+1}+T_{j-1}^{n+1}}{\Delta r^2}+\frac{T_{j+1}^{n+1}-T_{j-1}^{n+1}}{2r_j\Delta r}\right]
$$

The coefficient multiplying $T_{j-1}^{n+1}$ is

$$
-\frac{\alpha\Delta t}{2\Delta r^2} + \frac{\alpha\Delta t}{4r_j\Delta r}
$$

The coefficient multiplying $T_j^{n+1}$ is

$$
1+\frac{\alpha\Delta t}{\Delta r^2}
$$

The coefficient multiplying $T_{j+1}^{n+1}$ is

$$
-\frac{\alpha\Delta t}{2\Delta r^2} - \frac{\alpha\Delta t}{4r_j\Delta r}
$$

Therefore, the matrix coefficients are

$$
A_{j,j-1} = -\frac{\alpha\Delta t}{2\Delta r^2} + \frac{\alpha\Delta t}{4r_j\Delta r}
$$

$$
A_{j,j} = 1+\frac{\alpha\Delta t}{\Delta r^2}
$$

$$
A_{j,j+1} = -\frac{\alpha\Delta t}{2\Delta r^2} - \frac{\alpha\Delta t}{4r_j\Delta r}
$$

These correspond directly to the `lower`, `diagonal`, and `upper` coefficients constructed in `calculations.py`.

The right-hand side for an interior node is

$$
b_j = T_j^n + \frac{\alpha\Delta t}{2}\left[\frac{T_{j+1}^n-2T_j^n+T_{j-1}^n}{\Delta r^2}+\frac{T_{j+1}^n-T_{j-1}^n}{2r_j\Delta r}\right]+S\Delta t
$$

## 12. Crank-Nicolson Centreline Treatment

At the centreline, the cylindrical spatial operator becomes

$$
\left. \left( \frac{\partial^2 T}{\partial r^2} + \frac{1}{r}\frac{\partial T}{\partial r} \right) \right|_{r=0} \approx \frac{4(T_1-T_0)}{\Delta r^2}
$$

The governing equation at the centreline is therefore

$$
\frac{\partial T_0}{\partial t} = 4\alpha\frac{T_1-T_0}{\Delta r^2}+S
$$

Applying the Crank-Nicolson method gives

$$
\frac{T_0^{n+1}-T_0^n}{\Delta t} = \frac{1}{2}\left[4\alpha\frac{T_1^{n+1}-T_0^{n+1}}{\Delta r^2}+4\alpha\frac{T_1^n-T_0^n}{\Delta r^2}\right]+S
$$

Using

$$
Fo = \frac{\alpha\Delta t}{\Delta r^2}
$$

gives

$$
T_0^{n+1}-T_0^n = 2Fo\left[T_1^{n+1}-T_0^{n+1}+T_1^n-T_0^n\right]+S\Delta t
$$

Collecting the unknown terms at time level $n+1$ gives

$$
(1+2Fo)T_0^{n+1}-2FoT_1^{n+1} = T_0^n+2Fo(T_1^n-T_0^n)+S\Delta t
$$

Therefore, the first matrix row contains

$$
A_{0,0}=1+2Fo
$$

and

$$
A_{0,1}=-2Fo
$$

with the right-hand side

$$
b_0 = T_0^n+2Fo(T_1^n-T_0^n)+S\Delta t
$$

This allows the centreline temperature to evolve according to the governing heat equation while satisfying the cylindrical symmetry condition.

## 13. Crank-Nicolson Surface Boundary

The convective surface boundary condition is

$$
-k\left. \frac{\partial T}{\partial r} \right|_{r=R} = h(T_s-T_\infty)
$$

Using a backward difference,

$$
-k\frac{T_s^{n+1}-T_{N-2}^{n+1}}{\Delta r} = h(T_s^{n+1}-T_\infty)
$$

Expanding gives

$$
-\frac{k}{\Delta r}T_s^{n+1}+\frac{k}{\Delta r}T_{N-2}^{n+1} = hT_s^{n+1}-hT_\infty
$$

Rearranging,

$$
\frac{k}{\Delta r}T_{N-2}^{n+1}+\left(-\frac{k}{\Delta r}-h\right)T_s^{n+1} = -hT_\infty
$$

The final matrix row is therefore

$$
A_{N-1,N-2} = \frac{k}{\Delta r}
$$

$$
A_{N-1,N-1} = -\frac{k}{\Delta r}-h
$$

and

$$
b_{N-1} = -hT_\infty
$$

## 14. Matrix System

Combining the centreline treatment, interior-node equations, and surface boundary condition produces a linear system

$$
A\mathbf{T}^{n+1} = \mathbf{b}
$$

where

$$
\mathbf{T}^{n+1} = \begin{bmatrix} T_0^{n+1} \\ T_1^{n+1} \\ \vdots \\ T_{N-2}^{n+1} \\ T_{N-1}^{n+1} \end{bmatrix}
$$

The coefficient matrix has the general form

$$
A = \begin{bmatrix} 1+2Fo & -2Fo & 0 & \cdots & 0 \\ a_1 & b_1 & c_1 & \cdots & 0 \\ 0 & a_2 & b_2 & \ddots & \vdots \\ \vdots & \ddots & \ddots & \ddots & c_{N-2} \\ 0 & \cdots & 0 & \frac{k}{\Delta r} & -\frac{k}{\Delta r}-h \end{bmatrix}
$$

where the interior coefficients are

$$
a_j = -\frac{\alpha\Delta t}{2\Delta r^2}+\frac{\alpha\Delta t}{4r_j\Delta r}
$$

$$
b_j = 1+\frac{\alpha\Delta t}{\Delta r^2}
$$

$$
c_j = -\frac{\alpha\Delta t}{2\Delta r^2}-\frac{\alpha\Delta t}{4r_j\Delta r}
$$

The system is solved at every time step using

```python
T[n, :] = np.linalg.solve(A, b)