"""
This script simulates the cross section of two vortex rings next to each other 
on a 100x100 grid

Each vortex ring is represented by two vortices
The animation is saved as a gif file named vortex_animation.gif

The animation is not completely accurate and does not demonstrate leapfrogging behaiviour 
There is some math mistake in my code somewhere


@author: Javeria Rizwan
Feb 12, 2023
"""
import numpy as np
import matplotlib.pyplot as pl
import math
import pandas as pd
from matplotlib.animation import FuncAnimation

# function to convert cartesian coordinates to polar
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)


# lists to save xy velocities and the vortex center coordiantes for animation
vel_xs = []
vel_ys = []
x_vs = []
y_vs = []


dt = 2 ### timestep
Nsteps = 100 ### total number of timesteps

## Setting up initial conditions (vortex centres and circulation)

# Vortex rings
y_v = np.array([-40, 40, -50, 50]) ### insert the y-positions of the 4 vortices
x_v = np.array([-80, -80, -65, -65]) ### insert the x-positions of the 4 vortices
k = 1
k_v = np.array([-k, k, -k, k])  # sign of k corresponds to direction of flow
# magnitude of circulation 2pi*k has to be the same for all 


ngrid = 100 # dimension of simulation grid
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j] 

vel_x = np.zeros(np.shape(Y)) #this holds x-velocity
vel_y = np.zeros(np.shape(Y)) #this holds y-velocity

# masking radius for better visualization of the vortex centres
r_mask = 10 
#within this mask, you will not plot any streamline 
#so that you can see more clearly the movement of the vortex centres



# calculate initial velocity at each point at X and Y due to individual vortices
rows, columns = X.shape
for v in range(len(x_v)): #looping over each vortex
    print(f"Getting initial velocities due to vortex {v}")
    ### computing the total velocity field
    
    center = (x_v[v], y_v[v]) # center of this vortex
    k = k_v[v]
    
    
    ## convert the points on grid to the xy coordinates relative to vortex center
    xrel, yrel = X - center[0], Y - center[1] 
    
    r = np.sqrt(xrel**2 + yrel**2)
    
    ## calculate velocity contributions at each point on the grid due to this vortex
    ## This is the cartesian form of the velocity u = k/r in the phi direction
    ux = -k * (Y - center[1]) / r**2
    uy = k * (X - center[0]) / r**2
    
    vel_x += ux
    vel_y += uy
    
    ## lines for masking (set velocities to NaN)
    vel_x[r < r_mask] = np.nan
    vel_y[r < r_mask] = np.nan
    
    
# save velocity fields of initial conditions for animation
vel_xs.append(vel_x)
vel_ys.append(vel_y)
x_vs.append(x_v)
y_vs.append(y_v)


# Setting up the plot for initial condition
fig, ax = pl.subplots(1,1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, 'k+', markersize=10) 
# set up the boundaries of the simulation box
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines
ax.streamplot(X, Y, vel_x, vel_y, density=[1, 1]) 
pl.show()



# Evolution
count = 0
while count < Nsteps:
    
    ## Compute and update advection velocity from each vortex center
    adv_velocities_on_centers_x = np.zeros(x_v.shape)
    adv_velocities_on_centers_y = np.zeros(x_v.shape)
    
    x_v_new = x_v.copy()
    y_v_new = y_v.copy()
    
    
    for vort in range(len(x_v)):
        
        # center of this vortex
        i, j = x_v[vort], y_v[vort]
        
        # compute the total advection velocity on this center due to all the neigbouring vortices:
        for v in range(len(x_v)):
            center = (x_v[v], y_v[v])
            if vort == v: # the vortex does not contribute velocity to itself
                continue 
            k = k_v[v]
            
            # same math as above
            xrel, yrel = i - center[0], j - center[1] # cartesian coordinates of r vector
            r = np.sqrt(xrel**2 + yrel**2) # polar coordiante of r vector
            
            if r == 0:
                continue
            # velocity contribution due to this vortex
            u_y = -k * (j - center[1]) / r**2
            u_x = k * (i - center[0]) / r**2
            
            # update position of centers
            x_v_new[vort] += u_x*dt
            y_v_new[vort] += u_y*dt
            
    print("Old x", x_v, "New x", x_v_new)
    print("Old y", y_v, "New y", y_v_new)
    
    x_v, y_v = x_v_new, y_v_new
    
    # with new vortex centers, re-initialize the total velocity field
    vel_x = np.zeros(np.shape(Y)) #this holds x-velocity
    vel_y = np.zeros(np.shape(Y)) #this holds y-velocity
    
    for v in range(len(x_v)): #looping over each vortex
        print(f"Getting velocities at timestep {count} due to vortex {v}")
        
        center = (x_v[v], y_v[v]) # center of this vortex
        k = k_v[v]
        
        
        ## convert the points on grid to the xy coordinates relative to vortex center
        xrel, yrel = X - center[0], Y - center[1] 
        
        r = np.sqrt(xrel**2 + yrel**2)
        
        ## calculate velocity contributions at each point on the grid due to this vortex
        ux = -k * (Y - center[1]) / r**2
        uy = k * (X - center[0]) / r**2
        
        vel_x += ux
        vel_y += uy
        
        ## lines for masking (set velocities to NaN)
        vel_x[r < r_mask] = np.nan
        vel_y[r < r_mask] = np.nan
                
               
    # save velocity fields of conditions for animation
    vel_xs.append(vel_x)
    vel_ys.append(vel_y)
    x_vs.append(x_v)
    y_vs.append(y_v)
    
    count += 1


## animate the time evolution
fig = pl.figure()
ax = pl.axes(xlim=(-ngrid, ngrid), ylim=(-ngrid, ngrid))

def animate(i):
    # the following four lines clear out the previous streamlines
    for coll in ax.collections:
        coll.remove()
    for patch in ax.patches:
        patch.remove()
    
    vel_x, vel_y = vel_xs[i], vel_ys[i]
    x_v = x_vs[i]
    y_v = y_vs[i]
    
    ax.clear()
    ax.plot(x_v, y_v, 'k+', markersize=10) 
    ax.streamplot(X, Y, vel_x, vel_y, density=[1, 1])
    ax.set_title(f"Frame {i}")


anim = FuncAnimation(fig, animate, frames=len(vel_xs), interval=200, repeat=False)
anim.save('vortex_animation.gif', writer='pillow')
#anim.save('bar_ani.mp4', writer='ffmpeg')
#pl.show()