# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:43:27 2022

@author: ryust
"""

import numpy as np
import matplotlib.pyplot as plt
from bubble_ring_functions_ext import *

# Create 2 bubble rings
N     = 20  #vertices
nt    = 60  #time steps
theta = np.linspace(0,2*np.pi,N+1)
R     = 0.5   #radius
dt    = 0.01#time step
mind  = R * theta[1] #distance for resampling

#Bubble ring 1
r1_pos = np.array([R*np.cos(theta),R*np.sin(theta),np.zeros(N+1)]).T #first vertex repeat at end
r1_C   = 2*np.ones(N+1) #edges, repeat for first edge
r1_a   = 0.1*np.ones(N+1) #edges, repeat for first edge
r1_vol = volume(r1_pos,r1_a)
r1_N   = N

#Bubble ring 2
r2_pos = np.array([R*np.cos(theta) - 1.25,R*np.sin(theta) - 0.3, np.zeros(N+1)]).T #first vertex repeat at end
r2_C   = 2*np.ones(N+1) #edges, repeat for first edge
r2_a   = 0.1*np.ones(N+1) #edges, repeat for first edge
r2_vol = volume(r2_pos,r2_a)
r2_N   = N

#Plot to check initial conditions
# fig = plt.figure(figsize = (10,10))
# ax = plt.axes(projection='3d')
# plt.plot(r1_pos[:,0],r1_pos[:,1],r1_pos[:,2])
# ax.scatter(r1_pos[:,0],r1_pos[:,1],r1_pos[:,2],s = 100)

# plt.plot(r2_pos[:,0],r2_pos[:,1],r2_pos[:,2])
# ax.scatter(r2_pos[:,0],r2_pos[:,1],r2_pos[:,2],s = 100)

#%%
#Frames for plotting
frame   = [[np.copy(r1_pos), np.copy(r2_pos)]]
frame_a = [[np.copy(r1_a), np.copy(r2_a)]]

#Time integration
for t in range(nt):
    #Runge Kutta https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    # pos_new = pos_old + 1/6 * dt * (k1 + 2k2 + 2k3 + k4)
    # k1 = vel(t       ,pos_old)
    # k2 = vel(t + dt/2,pos_old + dt/2 * k1)
    # k3 = vel(t + dt/2,pos_old + dt/2 * k2)
    # k4 = vel(t + dt  ,pos_old + dt * k3)

    r1_posRK = np.copy(r1_pos)
    r2_posRK = np.copy(r2_pos)
    
    #Step 1
    r1_vel = velocity(r1_posRK,r1_C,r1_a,r2_pos,r2_C,r2_a)
    r1_K1  = dt * r1_vel
    r2_vel = velocity(r2_posRK,r2_C,r2_a,r1_pos,r1_C,r1_a)
    r2_K1  = dt * r2_vel
    
    #Step 2
    r1_pos_old       = np.copy(r1_posRK)
    r1_posRK[:-1,:] += 0.5 * r1_K1
    r1_posRK[-1,:]   = r1_posRK[0,:] #Last element is first point
    r2_pos_old       = np.copy(r2_posRK)
    r2_posRK[:-1,:] += 0.5 * r2_K1
    r2_posRK[-1,:]   = r2_posRK[0,:] #Last element is first point
    
    r1_a   = modify_thickness(r1_pos_old,r1_posRK,r1_a,r1_N) #Modify thickness over edges
    r2_a   = modify_thickness(r2_pos_old,r2_posRK,r2_a,r2_N) #Modify thickness over edges
    
    r1_vel = velocity(r1_posRK,r1_C,r1_a,r2_pos,r2_C,r2_a)
    r1_K2  = dt * r1_vel
    r2_vel = velocity(r2_posRK,r2_C,r2_a,r1_pos,r1_C,r1_a)
    r2_K2  = dt * r2_vel
        
    #Step 3
    r1_pos_old       = np.copy(r1_posRK)
    r1_posRK[:-1,:] += 0.5 * r1_K2
    r1_posRK[-1,:]   = r1_posRK[0,:] #Last element is first point
    r2_pos_old       = np.copy(r2_posRK)
    r2_posRK[:-1,:] += 0.5 * r2_K2
    r2_posRK[-1,:]   = r2_posRK[0,:] #Last element is first point
    
    r1_a   = modify_thickness(r1_pos_old,r1_posRK,r1_a,r1_N) #Modify thickness over edges
    r2_a   = modify_thickness(r2_pos_old,r2_posRK,r2_a,r2_N) #Modify thickness over edges
    
    r1_vel = velocity(r1_posRK,r1_C,r1_a,r2_pos,r2_C,r2_a)
    r1_K3  = dt * r1_vel
    r2_vel = velocity(r2_posRK,r2_C,r2_a,r1_pos,r1_C,r1_a)
    r2_K3  = dt * r2_vel
    
    #Step 4
    r1_pos_old       = np.copy(r1_posRK)
    r1_posRK[:-1,:] += r1_K3
    r1_posRK[-1,:]   = r1_posRK[0,:] #Last element is first point
    r2_pos_old       = np.copy(r2_posRK)
    r2_posRK[:-1,:] += r2_K3
    r2_posRK[-1,:]   = r2_posRK[0,:] #Last element is first point
    
    r1_a   = modify_thickness(r1_pos_old,r1_posRK,r1_a,r1_N) #Modify thickness over edges
    r2_a   = modify_thickness(r2_pos_old,r2_posRK,r2_a,r2_N) #Modify thickness over edges
    
    r1_vel = velocity(r1_posRK,r1_C,r1_a,r2_pos,r2_C,r2_a)
    r1_K4  = dt * r1_vel
    r2_vel = velocity(r2_posRK,r2_C,r2_a,r1_pos,r1_C,r1_a)
    r2_K4  = dt * r2_vel
    
    #RK4
    r1_pos[:-1,:] += (r1_K1 + 2*r1_K2 + 2*r1_K3 + r1_K4)/6
    r1_pos[-1,:]   = r1_pos[0,:]
    r2_pos[:-1,:] += (r2_K1 + 2*r2_K2 + 2*r2_K3 + r2_K4)/6
    r2_pos[-1,:]   = r2_pos[0,:]
    
    #Volume conservation (Inside air)
    r1_vol1 = volume(r1_pos,r1_a)
    r1_a    = r1_a * np.sqrt(r1_vol/r1_vol1)
    r2_vol1 = volume(r2_pos,r2_a)
    r2_a    = r2_a * np.sqrt(r2_vol/r2_vol1)
    
    #Resample
    r1_pos, r1_a, r1_C, r1_N = resample(r1_pos,r1_a,r1_C,mind)
    r2_pos, r2_a, r2_C, r2_N = resample(r2_pos,r2_a,r2_C,mind)

    #Burgers thickness flow
    r1_a[:-1] = burgers_flow(r1_pos,r1_C,r1_a,dt)
    r1_a[-1]  = r1_a[0]
    r2_a[:-1] = burgers_flow(r2_pos,r2_C,r2_a,dt)
    r2_a[-1]  = r2_a[0]
    
    #Remove sharp coners (hairpins)
    r1_pos, r1_a, r1_C, r1_N = hairpin_removal(r1_pos,r1_a,r1_C)
    r2_pos, r2_a, r2_C, r2_N = hairpin_removal(r2_pos,r2_a,r2_C)
    
    #Resample
    r1_pos, r1_a, r1_C, r1_N = resample(r1_pos,r1_a,r1_C,mind)
    r2_pos, r2_a, r2_C, r2_N = resample(r2_pos,r2_a,r2_C,mind)

    frame.append([np.copy(r1_pos),np.copy(r2_pos)])
    frame_a.append([np.copy(r1_a),np.copy(r2_a)])
    
#%%
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.set_box_aspect((2,2,1))

for t in range(nt+1):
    if t%15 == 0:
        ring1 = frame[t][0]
        plt.plot(ring1[:,0],ring1[:,1],ring1[:,2])
        ax.scatter(ring1[:,0],ring1[:,1],ring1[:,2],s = frame_a[t][0]*1000)
        
        ring2 = frame[t][1]
        plt.plot(ring2[:,0],ring2[:,1],ring2[:,2])
        ax.scatter(ring2[:,0],ring2[:,1],ring2[:,2],s = frame_a[t][1]*1000)