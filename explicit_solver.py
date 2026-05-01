import numpy as np 

import matplotlib.pyplot as pl 

import matplotlib.animation as animation 

 

# Physical parameters 

R = 0.001 #Radius of resistor (m) 

 

#Conduction term parameters 

k = 11.3 #thermal conductivity of resistor (W/mK) 

rho = 8400 #density 

c_p = 450 #specific heat capacity 

Q = 1000000 #heat generation per unit volume 

alpha = k / (rho*c_p) 

 

#Convection term parameters 

h = 50 #convection heat transfer coefficient (W/m^2K) 

T_infinity = 300 #ambient temperature (K) 

 

#Spatial node 

a = 0 

b = R 

dr = 0.0001 #distance between spatial nodes 

 

#Temporal node 

Duration = 130 #seconds 

dt = 0.001 #period between time nodes 

 

def PDEsolver(): 

 

#create spatial grid points 

 r = np.arange(0,b+dr,dr) 

 Nr = len(r) #number of nodes 

 

 #create temporal grid points 

 t = np.arange(0,Duration+dt, dt) 

 Nt = len(t) 

 

 #create the solution matrix 

 T = np.zeros((Nt,Nr)) 

 T[0,:] = T_infinity 

 

 #compute solution incrementally at each time step 

 for i in range(1,Nt): 

     #boundary conditions 

     #Symmetrical 

       T[i, 0] = T[i-1, 1] 

       T[i, -1] = (T[i-1, -2] + (h * dr / k) * T_infinity) / (1 + h * dr / k)  

      

     #interior nodes 

       for j in range(1, Nr-1): 

         term1 = (T[i-1, j+1] - 2*T[i-1, j] + T[i-1, j-1]) / dr**2 

         term2 = (T[i-1, j+1] - T[i-1, j-1]) / (2 * r[j] * dr) 

         source_term = Q / (rho * c_p) 

          

         T[i, j] = T[i-1, j] + dt * (alpha * (term1 + term2) + source_term) 

 

 return(T, r, t) 

 

T, r, t = PDEsolver() 

courant = alpha*(dt/dr**2) 

vmin = np.min(T) 

vmax = np.max(T) 

 

r_full = np.concatenate((-r[::-1][:-1], r))  # from -R to R (skip the duplicate zero) 

T_full = np.concatenate((T[-1, ::-1][:-1], T[-1, :]))  # mirror the temperature profile 

 

# Line graph 

#plot results 

pl.grid() 

pl.plot(r_full*1000, T_full) 

pl.xlabel("Radius (mm)") 

pl.ylabel("Temperature (K)") 

pl.title("Temperature Distribution") 

pl.show() 

 

print(courant) 

 

#Static plot 

#plotting polar heat map at the final time step 

T_final = T[-1, :] 

#define r and theta grids 

theta = np.linspace(0,2*np.pi,360) 

Theta_edges = np.linspace(0, 2*np.pi, len(theta)+1) 

R_edges = np.linspace(0, R, len(r)+1) 

R, Theta = np.meshgrid(R_edges, Theta_edges) 

#assigning the temperatures at the last time step to each r and theta points 

T_assign = np.tile(T_final, (len(Theta_edges)-1,1)) 

 

#convert polar coordinates to cartesian 

def polartocartesian(R,Theta): 

   x = R*np.cos(Theta) 

   y = R*np.sin(Theta) 

   return(x,y) 

 

x, y = polartocartesian(R,Theta) 

 

#Plot heatmap 

heatmap1 = pl.pcolormesh(x*1000, y*1000, T_assign, shading='auto', cmap='RdYlBu_r') 

pl.xlabel("x (mm)") 

pl.ylabel("y (mm)") 

pl.title("Temperature Distribution") 

pl.colorbar(heatmap1, label="Temperature (K)") 

#even axis 

pl.axis("equal") 

pl.show() 

  

#Animate the Temperature Distribution over the total duration 

#set up the first plot 

fig, ax = pl.subplots() 

heatmap2 = ax.pcolormesh(x*1000, y*1000, np.tile(T[0], (len(Theta_edges)-1,1)), shading='auto', cmap='RdYlBu_r', vmin=vmin, vmax=vmax) 

ax.set_xlabel("x (mm)") 

ax.set_ylabel("y (mm)") 

ax.set_title("Temperature Distribution") 

pl.colorbar(heatmap2, label="Temperature (K)") 

 

#repeat for consequent frames 

def animate(frames): 

   #replace current heatmap data with new data 

   T_assign = np.tile(T[frames], (len(Theta_edges)-1, 1)) 

   heatmap2.set_array(T_assign.ravel()) 

   ax.set_title(f"Temperature at t={t[frames]: .2f} s") 

   return[heatmap2] 

 

#Animation function 

heat_animation = animation.FuncAnimation(fig, animate, frames=range(0, len(t), 400), blit=False, repeat=False) 

heat_animation.save("Temperature_Distribution.mp4", writer='ffmpeg') 

 

pl.show() 
