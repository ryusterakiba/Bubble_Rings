# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:43:27 2022

@author: ryust
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.sparse as sp
import scipy.sparse.linalg as la

#%%
muRM = np.exp(-0.75)
dRM  = 0.5 * np.exp(0.25)
nu   = 1e-6
g    = np.array([0,0,-9.8])
At   = -1

# Biotsavart (run over points)
norm2 = lambda x: np.inner(x,x)
def biotsavartedge(p,v0,v1,C,a):
    r0 = v0 - p
    r1 = v1 - p
    T  = r1 - r0
    a2mu     = (a * muRM)**2
    crossr01 = np.cross(r0,r1)
    
    term1 = np.dot(r1,T) / (np.sqrt(a2mu + norm2(r1)) * (norm2(T)*a2mu + norm2(crossr01)))
    term2 = np.dot(r0,T) / (np.sqrt(a2mu + norm2(r0)) * (norm2(T)*a2mu + norm2(crossr01))) 
    return C / (4*np.pi) * (term1 - term2) * crossr01    
    
#Induction term (run over points)
def induction(v0,v1,v2,c0,c1,a0,a1):
    s0 = np.linalg.norm(v1 - v0)
    s1 = np.linalg.norm(v2 - v1)
    T0 = (v1 - v0)/s0
    T1 = (v2 - v1)/s1
    C  = (c0 + c1)/2
    
    kB      = 2 * np.cross(T0,T1) / (s0 + s1)
    logterm = np.log(s0 * s1 / (a0 * a1 * dRM**2))
    return C / (4*np.pi) * 0.5 * logterm * kB

# Boussinesq (run over edges, interpolate to points)
def boussinesq(v0,v1,C,a):
    edge = v1 - v0
    ds   = np.linalg.norm(edge)
    T    = edge / ds
    Atg_n= At * g - np.dot(At * g, T) * T
    
    coeff1 = 16 * np.pi**2 * nu * a**2 / (256 * np.pi**2 * nu**2 + C**2)
    coeff2 = np.pi * a**2 * C/ (256 * np.pi**2 * nu**2 + C**2)
    return coeff1 * Atg_n + coeff2 * np.cross(T,Atg_n)

def velocity(pos,N,C,a):
    #Calculates the velocity of the N vertices located at pos
    #Each edge has circulation C and thickness a
    pos_pad = np.vstack((pos[-2,:],pos)) #has one vertex before start
    point_vel = np.zeros((N,3))

    for p in range(N):
        #Over all vertices of filament
        for i in range(N): 
            #Biot Savart sum over edges (TERM 1 uBSdisc)
            point_vel[p,:] += biotsavartedge(pos[p,:],pos[i,:],pos[i+1,:],C[i],a[i])
        #Apply induction term  (TERM 2 uLIA)
        point_vel[p,:] += induction(pos_pad[p,:],pos_pad[p+1,:],pos_pad[p+2,:],C[i],C[i+1],a[p],a[p+1])
    
    #Boussinesq term (TERM 3) evaluated at edges then interpolated to vertices
    edge_vel = np.zeros((N,3))
    for i in range(N):
        #Boussinesq over edges
        edge_vel[i,:] = boussinesq(pos[i,:],pos[i+1,:],C[i],a[i])
    edge_vel = np.vstack((edge_vel[-1,:],edge_vel)) #pad one before start
    for p in range(N):
        #interpolate to vertices
        point_vel[p,:] += (edge_vel[p,:] + edge_vel[p+1,:]) / 2
    
    return point_vel

def modify_thickness(pos_old,pos_new,a,N):
    #Modify thickness over edges
    for p in range(N):
        l_old = np.linalg.norm(pos_old[p+1] - pos_old[p])
        l_new = np.linalg.norm(pos_new[p+1] - pos_new[p])  
        a[p] = a[p] * np.sqrt(l_old / l_new)
    a[-1] = a[0]
    return a

#Resampling
def resample(pos,a, N):
    # edge length and volume (cylinder)
    # interpolate curve
    edges = pos[1:] - pos[:-1]
    d = np.insert(np.cumsum([np.linalg.norm(edges[i,:]) for i in range(edges.shape[0])]),0,0)
    spaced = np.linspace(0,d.max(),N)
    
    #Interpolate
    pos_new = np.vstack((np.interp(spaced,d, pos[:,0]),np.interp(spaced,d, pos[:,1]),np.interp(spaced,d, pos[:,2]))).T
    a_new = np.interp(spaced,d,a)
    
    return pos_new,a_new

#Tangential flow (want to input all positions i guess)
def bergers_flow(pos,C,a,dt):
    edges = pos[1:] - pos[:-1]
    N     = edges.shape[0]
    ds    = np.array([np.linalg.norm(edges[i,:]) for i in range(N)])
    T     = np.array([edges[i,:] / ds[i] for i in range(N)])
    areas = np.pi * a**2
    
    #Compute Flux Eq 25 (phi * nu)
    flux  = np.zeros(N)
    grav  = np.array([np.dot(At * g, T[i]) for i in range(N)])
    coeff = 1 / (8 * np.pi)
    
    #some padding so indexing works out
    ds_pad   = np.hstack((ds[-1],ds))
    area_pad = np.hstack((areas[-1],areas))
    grav_pad = np.hstack((grav[-1],grav))
    ds_ave   = (ds_pad[1:] + ds_pad[:-1])/2
    
    for i in range(N): #Iterate over vertices
        pre = grav_pad[i] * area_pad[i]
        nex = grav_pad[i+1] * area_pad[i+1]
        
        if pre > max(0,-nex):
            flux[i] = coeff * grav_pad[i] * area_pad[i]**2
        elif nex < min(0,-pre):
            flux[i] = coeff * grav_pad[i+1] * area_pad[i+1]**2
        else:
            flux[i] = 0
        
    #Implicit Euler (Probably have to change to sparse matrices to keep numerical stability)    
    graph = np.diag(-np.ones(N),0) + np.diag(np.ones(N-1),1)
    graph[-1,0] = 1
    graph = sp.csr_matrix(graph)
    
    M    = sp.diags(ds)
    st   = sp.diags(ds_ave**-1 * C[:N]**2)
    L    = -graph.transpose().dot(st.dot(graph)) #basically -2 main diagonal, 1 off diagonal, 1 corners
    coef = 1 / (64*np.pi**2)
    nudt = nu / dt
    
    #Divided by dt, flux = phi * nu
    LHS  = nudt * M - 0.5 * coef * L
    RHS  = nudt * M.dot(areas[:N]) + graph.transpose().dot(flux)
    scl  = 1 / np.linalg.norm(RHS)
    # sol  = np.linalg.solve(LHS,RHS*scl)
    sol,info = la.cg(LHS,RHS*scl)
    if info != 0:
        print("oh no", info)
    sol  = sol/scl
    
    if any(sol < 0):
        print("negative encountered")
        return np.sqrt(np.abs(sol/np.pi))
    else:
        return np.sqrt(sol/np.pi)
    

#%%
# make vertices
N    = 20
nt   = 30
theta= np.linspace(0,2*np.pi,N+1)
z_pos= np.random.rand(N) * 0.0
z_pos= np.hstack((z_pos,z_pos[0]))
pos  = np.array([np.cos(theta),np.sin(theta),z_pos]).T #first vertex repeat at end
pos2 = np.vstack((pos[-2,:],pos)) #has one vertex before start
C    = np.ones(N+1) #edges, repeat for first edge
a    = 0.2*np.ones(N+1) #edges, repeat for first edge
dt   = 0.1

#%% 

frame = [np.copy(pos)]
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
    
    frame.append(np.copy(pos))
    
    #Bergers thickness flow
    a[:N] = bergers_flow(pos,C,a,dt)
    a[N]  = a[0]
    
#%%
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')

for t in range(nt+1):
    if t%10 == 0:
        plt.plot(frame[t][:,0],frame[t][:,1],frame[t][:,2])
