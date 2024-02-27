# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 12:54:25 2024

@author: yorub
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root

def eq(u, nu):
    #return nu**(3/2) * U**(3/2) * (1e6)*30 + (1e-3)*U**3*30*(173-(1e6*nu)/U) + 0.5*U**3 *30*15 - 23
    #return 1e4*np.sqrt(10)*nu*U**2 + 2e3*nu**(1/2)*U**(5/2) + 2e4**U**3 - 1e7*U**2*nu + 2.5e2*U**3 - 1.9e8
    #return 2*10**(5/2)*nu*U**2*(50+200) + 1/2 * U**3 * 50*200 - 4e8
    return 120*10**(5/2)*nu*u**2 + 250*u**3 - 2.5e8 

nu_water = 0.01
nus = np.linspace(1*nu_water, 1e4*nu_water, 501)

common_nus = {"Acetone": 0.00302,
              "Ethanol": 0.01074,
              "Mercury": 0.01526,
              "Glycerol": 14.12}

zeros = []

zeros_diff = []
for nu in nus:
    
    func = lambda U: eq(U, nu)
    #U_num = fsolve(func, 100)[0]
    U_num = root(func, 100).x[0]
    zeros.append(U_num)
    zeros_diff.append(100 - U_num )
    
plt.plot(nus, zeros)
plt.xlabel(r"$\nu$ (cm$^2$/s) ")
plt.ylabel("U (cm/s)")
colors = ["purple", "pink", "olive", "cyan"]

i=0
for liquid, nu in common_nus.items():
    plt.vlines(nu, 0, 100, label = liquid, colors=colors[i], linestyle = "dashed")
    i+=1
plt.legend()
plt.xscale("log")
plt.show()

plt.plot(nus, zeros_diff)
plt.xlabel(r"$\nu$ (cm$^2$/s) ")
plt.ylabel("Difference from U$_{in water}$ (%)")
plt.hlines(10, nus[0], nus[-1], label = "Professional", colors=colors[0], linestyle = "dashed")
plt.hlines(20,nus[0], nus[-1], label = "Recreational", colors=colors[1], linestyle = "dashed")
plt.legend()