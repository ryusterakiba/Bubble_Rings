# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 20:23:21 2022

@author: ryust
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la

#some constants
muRM = np.exp(-0.75)
dRM  = 0.5 * np.exp(0.25)
nu   = 1e-6
g    = np.array([0,0,-9.8])
At   = -1
hairpin_angle = -0.8
reconnection_length = 0.08
length_energy_factor = 20
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

def velocity(pos,C,a,pos_list,C_list,a_list):
    #Calculates the velocity of the N vertices located at pos
    #Each edge has circulation C and thickness a
    #pos_ext are positions of other bubble ring
    N       = pos.shape[0] - 1
    pos_pad = np.vstack((pos[-2,:],pos)) #has one vertex before start
    point_vel = np.zeros((N,3))

    for p in range(N):
        #Over all vertices of filament
        for f in range(len(pos_list)):
            #Over all vortex rings
            N_ext = pos_list[f].shape[0] - 1
            for i in range(N_ext):
                #Biot Savart sum over edges
                point_vel[p,:] += biotsavartedge(pos[p,:],pos_list[f][i,:],pos_list[f][i+1,:],
                                                  C_list[f][i],a_list[f][i])
            
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
    
#%% Bubble ring reconnection

def hairpin_removal(pos,a,C):
    #Remove sharp corners from a vortex filament
    N       = pos.shape[0] - 1    
    #Loop through pairs of edges and check for sharp angles (Hairpins)
    pos_new = [pos[-2]]
    a_new   = [a[-2]]
    change  = 0
    for j in range(N):
        e1 = pos[j] - pos_new[-1]
        if j == N-1:
            e2 = pos_new[1] - pos_new[0] 
        else:
            e2 = pos[j+1] - pos[j]
        de1= np.linalg.norm(e1)
        de2= np.linalg.norm(e2)
        cos_angle = np.dot(e1,e2) / de1 / de2
        if cos_angle < hairpin_angle: 
            #Remove hairpins by combining pair of edges
            #skip position
            #Adjust a
            change    = 1
            d_new     = np.linalg.norm(e1 + e2)
            if j == N-1:
                a_new[0] = 0.5 * (a_new[-1] + a_new[0]) * np.sqrt(d_new/(de1+de2))
                pos_new[0] = pos_new[-1]
            else:
                a_new[-1] = 0.5 * (a_new[-1] + a[j]) * np.sqrt(d_new/(de1+de2))     
        else:
            pos_new.append(pos[j])
            a_new.append(a[j])
    if change == 1:
        print('removed a hairpin')
        pos_new = np.array(pos_new)
        a_new   = np.array(a_new)
        C_new   = C[0] * np.ones(len(a_new))
        N_new   = pos_new.shape[0] - 1
    else:
        pos_new = pos
        a_new   = a
        C_new   = C
        N_new   = N
        
    return pos_new, a_new, C_new, N_new
            
# Reconnection helper functions
def evaluateIntegral(K,C,D,l):
    #Evaluate integral usign antiderivative
    L0  = C
    L1  = C+l
    SR0 = np.sqrt(K**2 + C**2 + D**2)
    SR1 = np.sqrt(K**2 + (C+l)**2 + D**2)
    lg = np.log(L1 + SR1) - np.log(L0 + SR0)

    if D > 1e-8:
        atan00 = np.arctan((K*L0)/(D*SR0))
        atan01 = np.arctan(L0/D)
        atan10 = np.arctan((K*L1)/(D*SR1))
        atan11 = np.arctan(L1/D)
        atn    = K*(atan10 - atan11 - atan00 + atan01)/D
        return lg + atn
    else:
        # D goes to 0
        f      = L0/(K + SR0) - L1/(K + SR1)
        return lg + f
     

def velocity_flux(eS, eT, p, a1, a2):
    # Evaluate velocity flux Apdx A1 Weissmann and Pinkall 2010 
    #eS, eT edge vectors, P is offset
    
    aa = a1*a2*muRM**2
    dS = np.linalg.norm(eS)
    dT = np.linalg.norm(eT)
    S  = eS / dS
    T  = eT / dT
    
    #Just in case?
    if (dS == 0 or dT == 0):
        return 0
    
    dp2  = np.linalg.norm(p)**2
    
    dSTp = np.linalg.norm(S + T)
    dSTm = np.linalg.norm(S - T)
    STp  = (S + T)/dSTp
    if dSTm == 0:
        STm  = (S - T)
    else:
        STm  = (S - T)/dSTm
    
    #H matrix
    pSTp = np.dot(p,STp)
    pSTm = np.dot(p,STm)
    fSTp = -pSTp/dSTp
    if dSTm == 0:
        fSTm = 0
    else:
        fSTm = -pSTm/dSTm
    H    = np.array([-pSTm,-pSTp,0])
    
    #Define K
    K2 = aa + dp2 - np.dot(H,H)
    K  = np.sqrt(max(0,K2))
    
    #G matrix
    G1 = np.array([dSTm/2, dSTp/2, 0])
    G2 = np.array([dSTm/2, -dSTp/2, 0])
    
    #Parallelogram
    A1 = np.array([-pSTm                            , -pSTp                             , 0])
    A2 = np.array([-pSTm + dS * dSTm/2              , -pSTp + dS * dSTp/2               , 0])
    A3 = np.array([-pSTm + dS * dSTm/2 + dT * dSTm/2, -pSTp + dS * dSTp/2 + dT * -dSTp/2, 0])
    A4 = np.array([-pSTm               + dT * dSTm/2, -pSTp               + dT * -dSTp/2, 0])
    
    #Parallel case
    MAX_ST = 1-1e-8
    ST  = max(-1,min(1,np.dot(S,T)))
    if (ST > MAX_ST or ST < -MAX_ST):
        #dS = L dT = l in Houdini code
        sgn = np.sign(ST)
        pT = np.dot(p,T)
        KK = aa + dp2 - pT**2
        
        XdTdS = pT + dT - sgn*dS
        XdT   = pT + dT
        XdS   = pT      - sgn*dS
        X0    = pT
        
        SRdTdS= np.sqrt(KK + XdTdS**2)
        SRdT  = np.sqrt(KK + XdT**2)
        SRdS  = np.sqrt(KK + XdS**2)
        SR0   = np.sqrt(KK + X0**2)

        fdTdS = SRdTdS - XdTdS * np.log(XdTdS + SRdTdS)
        fdT   = SRdT   - XdT   * np.log(XdT   + SRdT  )
        fdS   = SRdS   - XdS   * np.log(XdS   + SRdS  )
        f0    = SR0    - X0    * np.log(X0    + SR0   )
        
        parallel = sgn * (fdTdS - fdT - fdS + f0)
        return ST/(4*np.pi) * parallel
         
    C1 = np.dot(A1, G1)
    D1 = np.sqrt(max(0, np.dot(A1, A1) - C1*C1))
    
    C2 = np.dot(A2, G2)
    D2 = np.sqrt(max(0, np.dot(A2, A2) - C2*C2))
    
    C3 = -np.dot(A3, G1)
    D3 = np.sqrt(max(0, np.dot(A3, A3) - C3*C3))
    
    C4 = -np.dot(A4, G2)
    D4 = np.sqrt(max(0, np.dot(A4, A4) - C4*C4))
    
    #Compute Integral Eq 12
    I1 = evaluateIntegral(K, C1, D1, dS)
    I2 = evaluateIntegral(K, C2, D2, dT)
    I3 = evaluateIntegral(K, C3, D3, dS)
    I4 = evaluateIntegral(K, C4, D4, dT)

    I13= I1 - I3
    I24= I2 - I4
    I23= dS * I2 + dT * I3
    
    return ST/(4*np.pi) * (fSTp * (I24 + I13) + fSTm * (I24 - I13) + I23)

def quad_energy(p1, p2, p3, p4, a1, a2, a3, a4):
    #Kinetic Energy
    S1 = p2 - p1
    S2 = p3 - p2
    S3 = p4 - p3
    S4 = p1 - p4
    
    f11 = velocity_flux(S1, S1, p1 - p1, a1, a1)
    f12 = velocity_flux(S1, S2, p2 - p1, a1, a2)
    f13 = velocity_flux(S1, S3, p3 - p1, a1, a3)
    f14 = velocity_flux(S1, S4, p4 - p1, a1, a4)
    f1  = f11 + f12 + f13 + f14

    f21 = velocity_flux(S2, S1, p1 - p2, a2, a1)
    f22 = velocity_flux(S2, S2, p2 - p2, a2, a2)
    f23 = velocity_flux(S2, S3, p3 - p2, a2, a3)
    f24 = velocity_flux(S2, S4, p4 - p2, a2, a4)
    f2  = f21 + f22 + f23 + f24
    
    f31 = velocity_flux(S3, S1, p1 - p3, a3, a1)
    f32 = velocity_flux(S3, S2, p2 - p3, a3, a2)
    f33 = velocity_flux(S3, S3, p3 - p3, a3, a3)
    f34 = velocity_flux(S3, S4, p4 - p3, a3, a4)
    f3  = f31 + f32 + f33 + f34
    
    f41 = velocity_flux(S4, S1, p1 - p4, a4, a1)
    f42 = velocity_flux(S4, S2, p2 - p4, a4, a2)
    f43 = velocity_flux(S4, S3, p3 - p4, a4, a3)
    f44 = velocity_flux(S4, S4, p4 - p4, a4, a4)
    f4  = f41 + f42 + f43 + f44
    
    return 0.5 * np.sqrt(f1 + f2 + f3 + f4)

def neighbors(pos_list, a_list):
    #Find neighboring vertices 
    pos     = pos_list[0]
    a       = a_list[0]
    edges   = pos[1:] - pos[:-1]
    N       = edges.shape[0]
    ds      = np.array([np.linalg.norm(edges[i,:]) for i in range(N)])
    aveEdge = np.mean(ds)
    
    #Selection of near points parallel
    steps = np.floor( aveEdge / reconnection_length * 2.0)
    
    baseSearchDistance = reconnection_length
    if len(pos_list) == 1:
        #One bubble ring, check for splitting points
        nearpoints = [] #collect indices of points close together
        for i in range(N):
            #Iterate over all vertices
            nearpoints_i = [] #nearby points for vertex i
            searchDistance = max(baseSearchDistance, a[i]) + a[i]*1.5
            for j in range(N):
                distance = np.linalg.norm(pos[i,:] - pos[j,:])
                if distance < searchDistance + a[j]:
                    nearpoints_i.append(j)
            nearpoints.append(nearpoints_i)
    elif len(pos_list) == 2:
        #Two bubble rings, check for reconnection points
        other_pos  = pos_list[1]
        other_a    = a_list[1]
        nearpoints = [] #collect indices of points close together
        for i in range(N):
            #Iterate over all vertices
            nearpoints_i = [] #nearby points for vertex i
            searchDistance = max(baseSearchDistance, a[i]) + a[i]*1.5
            for j in range(other_pos.shape[0] - 1):
                distance = np.linalg.norm(pos[i,:] - other_pos[j,:])
                if distance < searchDistance + other_a[j]:
                    nearpoints_i.append(j)
            nearpoints.append(nearpoints_i)
    else:
        print('more than two rings')
        
    return nearpoints, steps
            
def is_reconnect(i1, i2, pos1, pos2, a1, a2, C1, C2, steps, nRings):
    if nRings == 1 and i1 == i2:
        # print('Same vertex')
        return 0, 0
    if (abs(C1[i1] - C2[i2]) > 1e-5 or abs(C1[i1]) < 1e-5 or abs(C2[i2]) < 1e-5):
        print('Different strength')
        return 0, 0
    
    N1 = pos1.shape[0] - 1
    N2 = pos2.shape[0] - 1 
    
    # print(np.linalg.norm(pos1[i1] - pos2[i2]))
    
    left_neighbors = np.mod(np.arange(i1 - steps,i1 + steps + 1), N1)
    right_neighbors= np.mod(np.arange(i2 - steps,i2 + steps + 1), N2)

    #Check if there is an overlap. not possible with 2 rings
    left_buffer = np.mod(np.arange(i1 - steps - 1 ,i1 + steps + 2), N1)
    right_buffer= np.mod(np.arange(i2 - steps - 1 ,i2 + steps + 2), N2)
    if nRings == 1:
        for i in range(len(left_buffer)):
            if any(right_buffer == left_buffer[i]):
                # print('Overlap')
                return 0, 0

    #Corner points
    PL1 = pos1[int(left_neighbors[0]),:]
    PL2 = pos1[int(left_neighbors[-1]),:]
    PR1 = pos2[int(right_neighbors[0]),:]
    PR2 = pos2[int(right_neighbors[-1]),:]
    
    Lnew = np.linalg.norm(PL1 - PR2) + np.linalg.norm(PL2 - PR1)
    Lold = np.linalg.norm(PL1 - PL2) + np.linalg.norm(PR1 - PR2)
    dL   = Lnew - Lold
    
    if dL > 0:
        # print('No length shortening')
        # print(dL)
        return 0, 0
    
    a12 = np.sqrt(a1[i1] * a2[i2])
    mC  = 0.5 * (C1[i1] + C2[i2])
    
    #Energy of neighboring quad
    dE = quad_energy(PL1, PL2, PR2, PR1, a1[i1], a12, a2[i2], a12)
    
    #Eq 10 of weissmann and pinkall 2010 for reconnection condition
    energy = length_energy_factor * dL + mC * dE
    
    if energy < 0:
        return 1, [left_neighbors, right_neighbors]
    # print('Energy high',energy)
    return 0, 0

def do_reconnect(left_neighbors, right_neighbors, pos1, pos2, a1, a2, C1, C2, nRings):
    #Do Splitting or Reconnection 
    if nRings == 1:
        #1 ring. pos1 = pos2, Split
        print('Splitting')
        pos = pos1
        a   = a1
        C   = C1
        N   = pos.shape[0] - 1
        
        i11 = int(left_neighbors[0])
        i12 = int(left_neighbors[-1])
        i21 = int(right_neighbors[0])
        i22 = int(right_neighbors[-1])
        
        #new positions
        if i21 > i12:
            idx1 = np.arange(i12,i21+1).astype(int)
        else:
            idx1 = np.concatenate((np.arange(i12,N), np.arange(0,i21+1))).astype(int)
        if i11 > i22:
            idx2 = np.arange(i22,i11+1).astype(int)
        else:
            idx2 = np.concatenate((np.arange(i22,N),np.arange(0,i11+1))).astype(int)
                                  
        posN1 = pos[idx1,:]
        posN2 = pos[idx2,:]
        
        aN1   = a[idx1]
        aN2   = a[idx2]
        #Average nearby a's for new edge
        aN1[-1] = 0.25 * (aN1[0] + aN1[-2] + aN1[-1] + a[int(left_neighbors[-2])])
        aN2[-1] = 0.25 * (aN2[0] + aN2[-2] + aN2[-1] + a[int(right_neighbors[-2])])
        
        #Shouldn't need to modify bc constant
        CN1   = C[idx1]
        CN2   = C[idx2]
        
        #Pad all with one element at end
        posN1 = np.vstack((posN1,posN1[0]))
        posN2 = np.vstack((posN2,posN2[0]))
        aN1   = np.hstack((aN1,aN1[0]))
        aN2   = np.hstack((aN2,aN2[0]))
        CN1   = np.hstack((CN1,CN1[0]))
        CN2   = np.hstack((CN2,CN2[0]))
        
        return [posN1,posN2], [aN1,aN2], [CN1,CN2]
        
    if nRings == 2:
        #2 rings, Reconnection
        print('Reconnecting')
        N1 = pos1.shape[0] - 1
        N2 = pos2.shape[0] - 1

        i11 = int(left_neighbors[0])
        i12 = int(left_neighbors[-1])
        i21 = int(right_neighbors[0])
        i22 = int(right_neighbors[-1])
        
        if i11 > i12:
            idx1 = np.arange(i12,i11+1).astype(int)
        else:
            idx1 = np.concatenate((np.arange(i12,N1), np.arange(0,i11+1))).astype(int)
        if i21 > i22:
            idx2 = np.arange(i22,i21+1).astype(int)
        else:
            idx2 = np.concatenate((np.arange(i22,N2), np.arange(0,i21+1))).astype(int)
            
        #Arange the new vertices
        pos = np.vstack((pos1[idx1],pos2[idx2],pos1[i12]))
        #Average nearby a's for new edge
        ai11 = 0.25 * (a1[i11] + a1[idx1[-2]] + a2[i22] + a2[int(right_neighbors[-2])])
        ai21 = 0.25 * (a2[i21] + a2[idx2[-2]] + a1[i12] + a1[int(left_neighbors[-2])])
        a   = np.hstack((a1[idx1[:-1]],ai11,a2[idx2[:-1]],ai21,a1[i12]))
        C   = np.hstack((C1[idx1],C2[idx2],C1[i12]))
        
        return [pos],[a],[C]
    else:
        print('cannot support more than 2 rings')
        return [pos1,pos2], [a1,a2], [C1,C2]


def Reconnection(pos_list,a_list,C_list):
    #We assume a max of two bubble rings at a time
    nRings = len(pos_list)
    
    #Get a list of nearby vertices from all primitives
    nearpoints, steps = neighbors(pos_list, a_list)
    # print(steps)
    
    #Check if reconnection 
    for i in range(len(nearpoints)):
        #For each edge of first bubble ring
        for j in range(len(nearpoints[i])):
            #For all the nearby points of that 
            isreconnect, neighbor = is_reconnect(i, nearpoints[i][j], pos_list[0], pos_list[nRings-1], 
                                                 a_list[0], a_list[nRings-1], C_list[0], C_list[nRings-1], steps, nRings)
            if isreconnect:
                pos_list,a_list,C_list = do_reconnect(neighbor[0], neighbor[1], pos_list[0], pos_list[nRings-1], 
                                                      a_list[0], a_list[nRings-1], C_list[0], C_list[nRings-1], nRings)
                
                return pos_list,a_list,C_list
        
    return pos_list,a_list,C_list
    
        
        
    





