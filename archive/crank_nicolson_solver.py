"""Problem Parameters""" 

# wire radius (m) 

wire_radius = 0.001  # roughly 6 gauge wire 

# Density (kg/m^3) 

rho = 8400 

# Specific heat capacity (j/kg) 

C = 450 

# Volumetric heat generation 

Q = 1e6 

S = Q/(rho*C)  # source term to be used in the equation 

# Initial Temperature (K) 

Ti = 300 

# Outside air temperaeture (K) 

Tinf = 300 

# Convection Coefficient 

h = 50  # taken for ambient air at 25C 

# Conduction Coefficient 

k = 11.3 

# thermal diffusitivity () 

alpha = k/(rho*C) 

 

 

def implicit_solver(r_steps, delta_t, t_end, error=0.01): 

    """Crank Nicolson Implementation of Parabolic heat diffusion PDE""" 

    delta_r = wire_radius/r_steps 

 

    t_steps = int(t_end/delta_t) 

    # results matrix 

    # [[r1,r2,r3], [r1,r2,r3],...] each row = new timestep 

    implicit_result = np.zeros((t_steps, r_steps)) 

    # create coeff matrix A 

    A = np.zeros((r_steps, r_steps)) 

    for index in range(0, r_steps): 

        # cyle through every row in coeff matrix 

        ri = index*delta_r  # +delta_r 

        # ai(Ti-1) + bi(Ti) + ci(Ti+1) = RHS 

        # due to polar coordinates, ai,bi,ci are functions of radius location so need to be recalculated for every location 

        if index != 0: 

            ai = ((alpha*delta_t)/(4*ri*delta_r)) - \ 

                ((alpha*delta_t)/(2*delta_r**2)) 

            bi = (1)+((alpha*delta_t)/(1*delta_r**2)) 

            ci = ((-alpha*delta_t)/(2*delta_r**2)) - \ 

                ((alpha*delta_t)/(4*ri*delta_r)) 

        if index == 0: 

            # Neumann B.C For Symmetry 

            A[index][0] = 1 

            A[index][1] = -1 

 

        elif index == r_steps-1: 

            # Robin B.C convection conduction 

            A[index][index] = -h-(k/delta_r) 

            A[index][index-1] = (k/delta_r) 

        else: 

            A[index][index-1] = ai 

            A[index][index] = bi 

            A[index][index+1] = ci 

    prev_T = np.ones(r_steps) 

    prev_T = prev_T * Ti 

    for time_step in (range(0, t_steps)): 

 

        b = np.zeros(r_steps) 

        for i in range(0, r_steps): 

            ri = i*delta_r 

            if i > 0 and i < r_steps-1: 

 

                b[i] = S*delta_t + prev_T[i] + ((alpha*delta_t)/2)(1/(2*ri*delta_r))( 

                    prev_T[i+1]-prev_T[i-1]) + ((alpha*delta_t)/2)(1/delta_r2)(prev_T[i+1]-2*prev_T[i]+prev_T[i-1]) 

 

        b[0] = 0  # Symmtery B.C 

        b[-1] = -h*Tinf  # Robin B.C 

 

        T = np.linalg.solve(A, b) 

        prev_T = np.copy(T) 

        implicit_result[time_step] = np.copy(T) 

    return implicit_result 

 

 

"""Implicit Solver Parameters""" 

# num r points 

r_steps = 50 

# num timesteps 

delta_t = 1e-3 

# end time (s) 

t_end = 10 

# matrix solve max error 

error = 0.01  # 1% error 

 

 

"""Get Solutions""" 

implicit_result = implicit_solver(r_steps, delta_t, t_end) 

# result is of the form [[T0,T1,T2...Tn-1,Tn],[]....] where each array is an individual timestep
