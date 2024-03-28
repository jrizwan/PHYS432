# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 20:42:37 2024

@author: yorub
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def advection_step(quantity, u, dt, dx):
    """
    Parameters
    ----------
    quantity : array
        density, momentum, or energy
    u : array
        velocity field.
    dt : float
        timestep.
    dx : float
        spatial step.

    Returns
    -------
    quantity : array
        density, momentum, or energy after one advection step.

    """
    
   
    J = np.zeros(n-1) # flux array
    
    for i in range(n-1):
        if u[i] > 0: # for positive velocity
            J[i] = quantity[i]*u[i] 
        else: # for negative velocity
            J[i] = quantity[i+1]*u[i] 
            
    # Update the quantity
    quantity[1:-1] = quantity[1:-1] - (dt/dx)*(J[1:]-J[:-1])
    quantity[0] = quantity[0] - (dt/dx)*J[0]
    quantity[-1] = quantity[-1] + (dt/dx)*J[-1]
    
    return quantity

gamma = 5/3 
n = 1000 # number of steps in simulation
dt = 1e-4

x = np.linspace(0, 1, n)  # Spatial grid
dx = x[1] - x[0]

### initial conditions
density = np.ones(n)  
momentum = np.zeros(n)
energy = np.ones(n)

P = (gamma - 1) * (energy - 0.5 * momentum ** 2 / density) ## pressure
speed_sound = gamma * P / density

# Gaussian perturbation in energy
energy = 1.0 + 40 * np.exp(-(x - 0.5) ** 2 / 0.0025)

#Mach number calculation
Mach  = momentum/density/np.sqrt(speed_sound)



#### Simulate #####
Densities, Machs = [], []
step = 0
while step < n:
    
    # compute advection velocity
    u = 0.5 * ((momentum[:-1] / density[:-1]) + (momentum[1:] / density[1:]))
    
    # Advect density and momentum
    density = advection_step(density, u, dt, dx)
    momentum = advection_step(momentum, u, dt, dx)
    
    # Compute pressure
    pressure = (gamma - 1) * (energy - 0.5 * momentum ** 2 / density)
    
    # evolve the momentum and calculate advection velocities
    momentum[1:-1] -= dt * (pressure[2:] - pressure[:-2]) / (2.0 * dx)
    momentum[0] -= dt * (pressure[1] - pressure[0]) / dx
    momentum[-1] -= dt * (pressure[-1] - pressure[-2]) / dx
    u = 0.5 * ((momentum[:-1] / density[:-1]) + (momentum[1:] / density[1:]))
    
    # Advect energy
    energy = advection_step(energy, u, dt, dx)
    
    # Compute pressure
    velocity = momentum / density
    pressure = (gamma - 1) * (energy - 0.5 * momentum ** 2 / density)
    
    # evlove energy energy
    energy[1:-1] -= dt * (pressure[2:] * velocity[2:] - pressure[:-2] * velocity[:-2]) / (2.0 * dx)
    energy[0] -= dt * (pressure[1] * velocity[1] - pressure[0] * velocity[0]) / dx
    energy[-1] -= dt * (pressure[-1] * velocity[-1] - pressure[-2] * velocity[-2]) / dx
    
    # calculate pressure and Mach number 
    pressure = (gamma - 1) * (energy - 0.5 * momentum ** 2 / density)
    speed_sound = np.sqrt(gamma * pressure / density)  # Sound speed
    Mach = velocity / speed_sound  
    
    
    Densities.append(np.copy(density))
    Machs.append(np.copy(Mach))
    step += 1

## animate the time evolution
fig, (ax1, ax2) = plt.subplots(2, 1)

def animate(i):
    # the following  lines clear out the previous streamlines
    for coll in ax1.collections:
        coll.remove()
    for patch in ax1.patches:
        patch.remove()
    for coll in ax2.collections:
        coll.remove()
    for patch in ax2.patches:
        patch.remove()
    
    density, Mach = Densities[i], Machs[i]
    
    ax1.clear()
    ax2.clear()
    ax1.plot(x, density)
    ax2.plot(x, Mach)
    
    ax1.set_ylim([0, 5])
    ax2.set_ylim([-5, 5])

    ax1.set_ylabel("Density")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Mach number")
    plt.suptitle(f"Frame {i}")
    


anim = FuncAnimation(fig, animate, frames=len(Densities), interval=75, repeat=False)
anim.save('DensitiesAndMachs.gif', writer='pillow')



#### Answering questions (see attached pdf on crowdmark)

## Find density ratio (Q1)
post_shock_max_density = np.max(density)  
density_ratio = post_shock_max_density / density[0]
print(f"Density ratio: {density_ratio}")

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(x, Densities[0])
ax2.plot(x, Densities[-1])
ax1.set_ylabel("Density pre-shock")
ax2.set_ylabel("Density post-shock")
ax2.set_xlabel("x")
plt.show()

## Estimate width of shock (Q2)
density, Mach = Densities[-1], Machs[-1]
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(x, density)
ax2.plot(x, Mach)

ax1.set_ylim([0, 5])
ax2.set_ylim([-5, 5])

ax1.set_ylabel("Density")
ax2.set_xlabel("x")
ax2.set_ylabel("Mach number")
plt.suptitle(f"After {n} steps")
plt.show()