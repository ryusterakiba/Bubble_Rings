# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:43:27 2022

@author: ryust
"""

import numpy as np
import matplotlib.pyplot as plt
from bubble_ring_functions_ext import *

# Create 2 bubble rings
N     = 40  #vertices
nt    = 30  #time steps
theta = np.linspace(0,2*np.pi,N+1)
R     = 0.5   #radius
dt    = 0.004#time step
mind  = R * theta[1] #distance for resampling

#Bubble ring 1
r1_pos = np.array([R*np.cos(theta),R*np.sin(theta),np.zeros(N+1)]).T #first vertex repeat at end
r1_C   = 2*np.ones(N+1) #edges, repeat for first edge
r1_a   = 0.04*np.ones(N+1) #edges, repeat for first edge
r1_vol = volume(r1_pos,r1_a)
r1_N   = N

#Bubble ring 2
N     = 40  #vertices
theta = np.linspace(0,2*np.pi,N+1)

r2_pos = np.array([R*np.cos(theta) - 1,R*np.sin(theta), np.zeros(N+1) - 0.2]).T #first vertex repeat at end
r2_C   = 2*np.ones(N+1) #edges, repeat for first edge
r2_a   = 0.04*np.ones(N+1) #edges, repeat for first edge
r2_vol = volume(r2_pos,r2_a)
r2_N   = N

pos_list = [r1_pos, r2_pos]
a_list   = [r1_a, r2_a]
C_list   = [r1_C,r2_C]
N_list   = [r1_N,r2_N]
vol_list = [r1_vol,r2_vol]

#%% Plot to check initial conditions
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.set_box_aspect(aspect = (2,1,1))

plt.plot(r1_pos[:,0],r1_pos[:,1],r1_pos[:,2])
ax.scatter(r1_pos[:,0],r1_pos[:,1],r1_pos[:,2],s = 100)

plt.plot(r2_pos[:,0],r2_pos[:,1],r2_pos[:,2])
ax.scatter(r2_pos[:,0],r2_pos[:,1],r2_pos[:,2],s = 100)
#%%
#Frames for plotting
frame   = [np.copy(pos_list)]
frame_a = [np.copy(a_list)]

#Time integration
for t in range(nt):
    #Runge Kutta https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    # pos_new = pos_old + 1/6 * dt * (k1 + 2k2 + 2k3 + k4)
    # k1 = vel(t       ,pos_old)
    # k2 = vel(t + dt/2,pos_old + dt/2 * k1)
    # k3 = vel(t + dt/2,pos_old + dt/2 * k2)
    # k4 = vel(t + dt  ,pos_old + dt * k3)

    nRings = len(pos_list)
    posRK  = np.copy(pos_list)

    K1 = []
    K2 = []
    K3 = []
    K4 = []
    for i in range(nRings):
        #Step 1
        vel = velocity(posRK[i],C_list[i],a_list[i],pos_list,C_list,a_list)
        K1.append(dt * vel)
       
        #Step 2
        pos_old          = np.copy(posRK[i])
        posRK[i][:-1,:] += 0.5 * K1[i]
        posRK[i][-1,:]   = posRK[i][0,:] #Last element is first point

        a_list[i] = modify_thickness(pos_old,posRK[i],a_list[i],N_list[i]) #Modify thickness over edges

        vel = velocity(posRK[i],C_list[i],a_list[i],pos_list,C_list,a_list)
        K2.append(dt * vel)
  
        #Step 3
        pos_old          = np.copy(posRK[i])
        posRK[i][:-1,:] += 0.5 * K2[i]
        posRK[i][-1,:]   = posRK[i][0,:] #Last element is first point
    
        a_list[i] = modify_thickness(pos_old,posRK[i],a_list[i],N_list[i]) #Modify thickness over edges

        vel = velocity(posRK[i],C_list[i],a_list[i],pos_list,C_list,a_list)
        K3.append(dt * vel)

        #Step 4
        pos_old          = np.copy(posRK[i])
        posRK[i][:-1,:] += K3[i]
        posRK[i][-1,:]   = posRK[i][0,:] #Last element is first point
    
        a_list[i] = modify_thickness(pos_old,posRK[i],a_list[i],N_list[i]) #Modify thickness over edges
        
        vel = velocity(posRK[i],C_list[i],a_list[i],pos_list,C_list,a_list)
        K4.append(dt * vel)
    
    
    for i in range(nRings):
    #RK4
        pos_list[i][:-1,:] += (K1[i] + 2*K2[i] + 2*K3[i] + K4[i])/6
        pos_list[i][-1,:]   = pos_list[i][0,:]

    #Volume conservation (Inside air)

        vol = volume(pos_list[i],a_list[i])
        a_list[i] = a_list[i] * np.sqrt(vol_list[i]/vol)
    
    #Resample
        pos_list[i], a_list[i], C_list[i], N_list[i] = resample(pos_list[i],a_list[i],C_list[i],mind)

    #Burgers thickness flow
        a_list[i][:-1] = burgers_flow(pos_list[i],C_list[i], a_list[i],dt)
        a_list[i][-1]  = a_list[i][0]

    
    # Hairpin and Reconnection
    for i in range(nRings):
    #Remove sharp coners (hairpins)
        pos_list[i], a_list[i], C_list[i], N_list[i] = hairpin_removal(pos_list[i],a_list[i],C_list[i])
        
    # Reconnection
    pos_list,a_list,C_list = Reconnection(pos_list,a_list,C_list)
    nRings = len(pos_list)
    
    # Resample
    for i in range(nRings):
        pos_list[i], a_list[i], C_list[i], N_list[i] = resample(pos_list[i],a_list[i],C_list[i],mind)
        

    frame.append(np.copy(pos_list))
    frame_a.append(np.copy(a_list))

    
#%%
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.set_box_aspect(aspect = (2,1,1))
# ax.view_init(170, 60)
ax.view_init(-160, 100)
ax.invert_zaxis()    

for t in range(nt+1):
# for t in range(14,25):
    if t%5 == 0:
        for i in range(len(frame[t])):
            ring = frame[t][i]
            plt.plot(ring[:,0],ring[:,1],ring[:,2])
            ax.scatter(ring[:,0],ring[:,1],ring[:,2],s = frame_a[t][i][0]*1000)











