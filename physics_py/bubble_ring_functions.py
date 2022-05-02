# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 20:23:21 2022

@author: ryust
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la

muRM = np.exp(-0.75)
dRM  = 0.5 * np.exp(0.25)
nu   = 1e-6
g    = np.array([0,0,-9.8])
At   = -1

#%% Velocity solver

def biotsavartedge(p,v0,v1,C,a):
    # Biotsavart (run over points)
    norm2    = lambda x: np.inner(x,x)
    r0       = v0 - p
    r1       = v1 - p
    T        = r1 - r0
    a2mu     = (a * muRM)**2
    crossr01 = np.cross(r0,r1)
    
    term1 = np.dot(r1,T) / (np.sqrt(a2mu + norm2(r1)) * (norm2(T)*a2mu + norm2(crossr01)))
    term2 = np.dot(r0,T) / (np.sqrt(a2mu + norm2(r0)) * (norm2(T)*a2mu + norm2(crossr01))) 
    return C / (4*np.pi) * (term1 - term2) * crossr01    
    
def induction(v0,v1,v2,c0,c1,a0,a1):
    #Induction term (run over points)
    s0 = np.linalg.norm(v1 - v0)
    s1 = np.linalg.norm(v2 - v1)
    T0 = (v1 - v0)/s0
    T1 = (v2 - v1)/s1
    C  = (c0 + c1)/2
    
    kB      = 2 * np.cross(T0,T1) / (s0 + s1)
    logterm = np.log(s0 * s1 / (a0 * a1 * dRM**2))
    return C / (4*np.pi) * 0.5 * logterm * kB

def boussinesq(v0,v1,C,a):
    # Boussinesq (run over edges, interpolate to points)
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
        point_vel[p,:] += induction(pos_pad[p,:],pos_pad[p+1,:],pos_pad[p+2,:],C[p],C[p+1],a[p],a[p+1])
    
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

#%% Functions for stabalizing integration
def modify_thickness(pos_old,pos_new,a,N):
    #Modify thickness over edges
    for p in range(N):
        l_old = np.linalg.norm(pos_old[p+1] - pos_old[p])
        l_new = np.linalg.norm(pos_new[p+1] - pos_new[p])
        a[p] = a[p] * np.sqrt(l_old / l_new)
    a[-1] = a[0]
    return a

def resample(pos,a,C,mind):
    # Resampling
    # for position, thickness, circulation
    # to keep distance of mind between points
    edges = pos[1:] - pos[:-1]
    N     = edges.shape[0]
    d     = np.insert(np.cumsum([np.linalg.norm(edges[i,:]) for i in range(N)]),0,0)
    newN  = int(d.max() // mind)
    spaced= np.linspace(0,d.max(),newN)
    
    #Interpolate from points at location "d" to points at "spaced"
    pos_new = np.vstack((np.interp(spaced,d, pos[:,0]),np.interp(spaced,d, pos[:,1]),np.interp(spaced,d, pos[:,2]))).T
    a_new   = np.sqrt(np.interp(spaced,d,a**2))
    C_new   = np.interp(spaced,d,C)
    
    return pos_new,a_new,C_new,newN - 1

def volume(pos,a):
    #Volume of bubble ring 
    # keep constant with initial by scaling thickness
    edges = pos[1:] - pos[:-1]
    ds    = np.array([np.linalg.norm(edges[i,:]) for i in range(edges.shape[0])])
    vol   = np.dot(ds,a[:-1]**2) * np.pi
    return vol

#%% Tangential flow component 
def burgers_flow(pos,C,a,dt):
    #Tangential flow (thickness update)
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
        
    #Implicit Euler (sparse matrices to keep numerical stability)    
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
    
