# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:43:27 2022

@author: ryust
"""

import numpy as np
import matplotlib.pyplot as plt
from bubble_ring_functions import *

#%%
# make vertices
N    = 20
nt   = 40
R    = 0.5
theta= np.linspace(0,2*np.pi,N+1)
z_pos= np.zeros(N+1)
pos  = np.array([R*np.cos(theta),R*np.sin(theta),z_pos]).T #first vertex repeat at end
C    = 2*np.ones(N+1) #edges, repeat for first edge
a    = 0.1*np.ones(N+1) #edges, repeat for first edge
dt   = 0.05
mind = R * theta[1]
vol0 = volume(pos,a)

#Tilt the ring
tilt = -10 * np.pi/180
rot  = np.array([[1, 0, 0],
                 [0, np.cos(tilt), -np.sin(tilt)],
                 [0, np.sin(tilt), np.cos(tilt)]])
pos_tilt = np.dot(rot,pos.T)
pos_tilt = pos_tilt.T
pos = pos_tilt

# fig = plt.figure(figsize = (10,10))
# ax = plt.axes(projection='3d')
# plt.plot(pos[:,0],pos[:,1],pos[:,2])
# ax.scatter(pos[:,0],pos[:,1],pos[:,2],s = 100)

# plt.plot(pos_tilt[:,0],pos_tilt[:,1],pos_tilt[:,2])
# ax.scatter(pos_tilt[:,0],pos_tilt[:,1],pos_tilt[:,2],s = 100)
#%% 
frame = [np.copy(pos)]
frame_a = [np.copy(a)]
for t in range(nt):
    #Runge Kutta https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    # pos_new = pos_old + 1/6 * dt * (k1 + 2k2 + 2k3 + k4)
    # k1 = vel(t       ,pos_old)
    # k2 = vel(t + dt/2,pos_old + dt/2 * k1)
    # k3 = vel(t + dt/2,pos_old + dt/2 * k2)
    # k4 = vel(t + dt  ,pos_old + dt * k3)

    pos_RK = np.copy(pos)
    
    #Step 1
    point_vel = velocity(pos_RK,N,C,a)
    K1        = dt * point_vel
    
    #Step 2
    pos_old       = np.copy(pos_RK)
    pos_RK[:N,:] += 0.5 * K1
    pos_RK[N,:]   = pos_RK[0,:] #Last element is first point
    
    a = modify_thickness(pos_old,pos_RK,a,N) #Modify thickness over edges
    point_vel = velocity(pos_RK,N,C,a)
    K2        = dt * point_vel
        
    #Step 3
    pos_old       = np.copy(pos_RK)
    pos_RK[:N,:] += 0.5 * K2
    pos_RK[N,:]   = pos_RK[0,:] #Last element is first point
    
    a = modify_thickness(pos_old,pos_RK,a,N) #Modify thickness over edges
    point_vel = velocity(pos_RK,N,C,a)
    K3        = dt * point_vel
    
    #Step 4
    pos_old       = np.copy(pos_RK)
    pos_RK[:N,:] += K3
    pos_RK[N,:]   = pos_RK[0,:] #Last element is first point
    
    a = modify_thickness(pos_old,pos_RK,a,N) #Modify thickness over edges
    point_vel = velocity(pos_RK,N,C,a)
    K4        = dt * point_vel
    
    #RK4
    pos[:N,:] += (K1 + 2*K2 + 2*K3 + K4)/6
    pos[N,:]   = pos[0,:]
    
    #Volume conservation (Inside air)
    vol1 = volume(pos,a)
    a    = a * np.sqrt(vol0/vol1)
    
    #Resample
    pos, a, C, N = resample(pos,a,C,mind)

    #Burgers thickness flow
    a[:N] = burgers_flow(pos,C,a,dt)
    a[N]  = a[0]

    frame.append(np.copy(pos))
    frame_a.append(np.copy(a))
    
#%%
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')

for t in range(nt+1):
    if t%10 == 0:
        plt.plot(frame[t][:,0],frame[t][:,1],frame[t][:,2])
        ax.scatter(frame[t][:,0],frame[t][:,1],frame[t][:,2],s = frame_a[t]*1000)
        
        
# plt.savefig("bubble_ring_py.png",dpi=300)
