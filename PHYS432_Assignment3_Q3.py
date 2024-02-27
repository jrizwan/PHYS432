"""
This script simulates the flow of lava on a one-dimensional grid
The lava flows on an inclined plane with angle alpha

This assumes viscocity of 10^3 cm^2 / second as estimated in class

@author: Javeria Rizwan
Feb 21, 2023
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from matplotlib.animation import FuncAnimation
import scipy

alpha = np.pi/6 # the angle of the inclined plane
g = 9.8
const = g*np.sin(alpha)

## parameters of analytic solution
Ngrid = 100
Nsteps = 150
dt = 5e-4
dx = 5e-2

V = 10**3 # viscocity estimated in class

beta = V*dt/dx**2

## Tri-diagonal Matrix for analytic solution
n=Ngrid
A = np.eye(n) * (1.0 + 2.0 * beta) + np.eye(n, k=1)*-beta + np.eye(n, k=-1)*-beta


us = [] # save the velocity grid at each timestep 



x = np.arange(0, Ngrid*dx, dx) # Spatial grid
v = np.ones(x.shape)*const*dt*-1 # vector with constant offset for gravity

# u is velocity at each spatial point
u = np.zeros(x.shape) # start fluid at rest everywhere 
us.append(np.copy(u))


## Stress-free condition
A[-1][-1] = 1 + beta

# evovlve through time
for ct in range(Nsteps):
    T = u - v
    u = scipy.linalg.solve(A, T)
    ### No slip boundary condition: 
    u[0] = 0
    
    us.append(np.copy(u))

    

## analytic steady state ##
H = x[-1] # the edge of the lava
u_analytic = -const/V * (1/2 * x**2 - H*x)


## animate the time evolution
fig = plt.figure()
ax = plt.axes()

def animate(i):
    # the following four lines clear out the previous streamlines
    for coll in ax.collections:
        coll.remove()
    for patch in ax.patches:
        patch.remove()
    
    u = us[i]
    
    
    ax.clear()
    ax.plot(x, u, label = "numerical")
    ax.plot(x, u_analytic, label = "steady state analytic")
    
    ax.set_xlabel("x")
    ax.set_ylabel("Velocity")
    ax.set_title(f"Frame {i}")
    ax.legend()


anim = FuncAnimation(fig, animate, frames=len(us), interval=200, repeat=False)
anim.save('Lava_WithAnalyticComparison.gif', writer='pillow')
    