#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 17:53:15 2020

@author: adrianzebrowski
"""
import numpy as np 
    
def maccormackFD(case,c_max,dx,gamma):
    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
        t_final = 0.035
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
    
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/2))
    x = np.linspace(0,L,ix)
    
    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_pred = np.zeros((3,ix))
    Q_update = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_pred = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        for i in range(1,ix-1): # PREDICTOR
            Q_pred[:,i] = Q[:,i]-(dt/dx)*(F[:,i+1]-F[:,i]) # calculate predictor value of Q
            
            rho_pred = Q_pred[0,i]
            u_pred = Q_pred[1,i]/rho_pred
            e_pred = Q_pred[2,i]/rho_pred-u_pred**2/2
            p_pred = e_pred*(gamma-1)*rho_pred # calculate values of rho, u, e, and p based on predictor Q
            
            F_pred[0,i] = rho_pred*u_pred # calculate predictor flux
            F_pred[1,i] = rho_pred*u_pred**2+p_pred
            F_pred[2,i] = rho_pred*u_pred*(e_pred+p_pred/rho_pred+u_pred**2/2)
            
        F_pred[:,0] = F[:,0]
        F_pred[:,ix-1] = F[:,ix-1]
        
        for i in range(1,ix-1): # CORRECTOR
            Q_update[:,i] = 0.5*(Q[:,i]+Q_pred[:,i]-(dt/dx)*(F_pred[:,i]-F_pred[:,i-1]))
        
        Q[:,1:ix-1] = Q_update[:,1:ix-1]
        rho = Q[0,:]
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
    
        F[0,:] = rho*u
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        
    return rho, u, p, e, x

def lax_fredrich(case,c_max,dx,gamma):
    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
        t_final = 0.035
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
    
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/2))
    x = np.linspace(0,L,ix)
    
    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_update = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    s = dx/dt
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        for i in range(1,ix-1): # F_half    
            F_half[:,i] = 0.5*(F[:,i]+F[:,i+1]-s*(Q[:,i+1]-Q[:,i])) # calculate F_half
        
        F_half[:,0] = F[:,0]
        F_half[:,ix-1] = F[:,ix-1]
        
        for j in range(1,ix-1): # Q_update
            Q_update[:,j] = Q[:,j]-(dt/dx)*(F_half[:,j]-F_half[:,j-1])
        
        Q[:,1:ix-1] = Q_update[:,1:ix-1]
        rho = Q[0,:]
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
    
        F[0,:] = rho*u
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        s = dx/dt
        
    return rho, u, p, e, x

def rusanov(case,c_max,dx,gamma):
    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
        t_final = 0.035
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
    
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/2))
    x = np.linspace(0,L,ix)
    
    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_update = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        for i in range(0,ix): # F_half  
                
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
                    
            s = np.maximum(np.abs(u[i])+a[i],np.abs(u[ip])+a[ip])
            F_half[:,i] = 0.5*(F[:,i]+F[:,ip]-s*(Q[:,ip]-Q[:,i])) # calculate F_half
        
        for i in range(0,ix): # Q_update
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
                    
            Q_update[:,i] = Q[:,i]-(dt/dx)*(F_half[:,i]-F_half[:,im])
        
        Q = Q_update[:,:]
        rho = Q[0,:]
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
    
        F[0,:] = rho*u
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        s = dx/dt
        
    return rho, u, p, e, x

def maccormackFV(case,c_max,dx,gamma):
    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
        t_final = 0.035
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
    
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/2))
    x = np.linspace(0,L,ix)
    
    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_star = np.zeros((3,ix))
    Q_starstar = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    F_halfstar = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        for i in range(1,ix-1): # F_half calculation
            F_half[:,i] = F[:,i+1]
        
        F_half[:,0] = F[:,0] # apply boundary conditions for F_half
        F_half[:,ix-1] = F[:,ix-1]
        
        for i in range(1,ix-1): # Q_star
            Q_star[:,i] = Q[:,i]-(dt/dx)*(F_half[:,i]-F_half[:,i-1]) 
            
        Q_star[:,0] = Q[:,0] # Apply boundary conditions for Q_star
        Q_star[:,ix-1] = Q[:,ix-1]
        
        rho_star = Q_star[0,:] # Calculate star variables
        u_star = Q_star[1,:]/rho_star
        e_star = Q_star[2,:]/rho_star-u_star**2/2
        p_star = e_star*(gamma-1)*rho_star
    
        F_halfstar[0,:] = rho_star*u_star # Calculate star values of F
        F_halfstar[1,:] = rho_star*u_star**2+p_star
        F_halfstar[2,:] = rho_star*u_star*(e_star+p_star/rho_star+u_star**2/2)
        
        for i in range(1,ix-1): # Q_starstar
            Q_starstar[:,i] = Q[:,i]-(dt/dx)*(F_halfstar[:,i]-F_halfstar[:,i-1])
            
        Q_starstar[:,0] = Q[:,0] # Apply boundary conditions for Q_starstar
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        Q = 0.5*(Q_star+Q_starstar) # Average Qstar and Qstarstar for next timelevel value of Q
        
        rho = Q[0,:] # Calculate all variables from Q
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
        
        F[0,:] = rho*u # Calculate new flux vector 
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a) # Advance timestep
        t = t+dt
        
    return rho, u, p, e, x

def maccormackFCT(case,c_max,dx,gamma):
    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
        t_final = 0.035
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
    
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/2))
    x = np.linspace(0,L,ix)
    
    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_star = np.zeros((3,ix))
    Q_starstar = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_halfstarL = np.zeros((3,ix))
    F_halfstarH = np.zeros((3,ix))
    F_halfL = np.zeros((3,ix))
    F_halfH = np.zeros((3,ix))
    A_half = np.zeros((3,ix))
    A_halfc = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        ### PREDICTOR STEP
        for i in range(1,ix-1): # F_halfH calculation using higher order MacCormack method
            F_halfH[:,i] = F[:,i+1]
            
        F_halfH[:,0] = F[:,0] # apply boundary conditions for F_half
        F_halfH[:,ix-1] = F[:,ix-1]
            
        for i in range(1,ix-1): # F_halfL calculation using lower order Rusanov method 
            s = np.maximum(np.abs(u[i])+a[i],np.abs(u[i+1])+a[i+1])
            F_halfL[:,i] = 0.5*(F[:,i]+F[:,i+1]-s*(Q[:,i+1]-Q[:,i])) # calculate F_half
            
        F_halfL[:,0] = F[:,0]
        F_halfL[:,ix-1] = F[:,ix-1]
            
        A_half = F_halfH-F_halfL # antidiffusion flux
            
        for i in range(1,ix-1): # intermediate MP solution using only low order F
            Q_star[:,i] = Q[:,i]-(dt/dx)*(F_halfL[:,i]-F_halfL[:,i-1]) 
                
        Q_star[:,0] = Q[:,0] # Apply boundary conditions for Q_star
        Q_star[:,ix-1] = Q[:,ix-1]
            
        for i in range(1,ix-2): # Calculate corrected antidiffusion flux
            A_halfc[0,i] = np.sign(A_half[0,i])*max(0,min(np.abs(A_half[0,i]),np.sign(A_half[0,i])*(Q_star[0,i+2]-Q_star[0,i+1])*dx/dt,np.sign(A_half[0,i])*(Q_star[0,i]-Q_star[0,i-1])*dx/dt))
            A_halfc[1,i] = np.sign(A_half[1,i])*max(0,min(np.abs(A_half[1,i]),np.sign(A_half[1,i])*(Q_star[1,i+2]-Q_star[1,i+1])*dx/dt,np.sign(A_half[1,i])*(Q_star[1,i]-Q_star[1,i-1])*dx/dt))
            A_halfc[2,i] = np.sign(A_half[2,i])*max(0,min(np.abs(A_half[2,i]),np.sign(A_half[2,i])*(Q_star[2,i+2]-Q_star[2,i+1])*dx/dt,np.sign(A_half[2,i])*(Q_star[2,i]-Q_star[2,i-1])*dx/dt))
            
        A_halfc[:,0] = A_half[:,0] # Need to apply boundary conditions on the corrected antidiffusion flux
        A_halfc[:,ix-2] = A_half[:,ix-2]
        A_halfc[:,ix-1] = A_half[:,ix-1]
    
        for i in range(1,ix-1): # flux corrected solution for predictor
            Q_star[:,i] =  Q_star[:,i]-(dt/dx)*(A_halfc[:,i]-A_halfc[:,i-1]) 
       
        Q_star[:,0] = Q[:,0] # Apply boundary conditions for Q_star again just to be safe while I debug this
        Q_star[:,ix-1] = Q[:,ix-1]
        
        rho_star = Q_star[0,:] # Calculate star variables from the corrected Q_star predictor
        u_star = Q_star[1,:]/rho_star
        e_star = Q_star[2,:]/rho_star-u_star**2/2
        p_star = e_star*(gamma-1)*rho_star
        a_star = np.sqrt(gamma*p_star/rho_star)
        
        ### END OF PREDICTOR STEP
        
        ### CORRECTOR STEP
    
        F_halfstarH[0,:] = rho_star*u_star # Calculate star values of Fhalf using higher order MacCormack method
        F_halfstarH[1,:] = rho_star*u_star**2+p_star
        F_halfstarH[2,:] = rho_star*u_star*(e_star+p_star/rho_star+u_star**2/2)
        
        for i in range(1,ix-1): # F_halfstarL calculation using lower order Rusanov method 
            s = np.maximum(np.abs(u[i])+a[i],np.abs(u[i+1])+a[i+1])
            F_halfstarL[:,i] = 0.5*(F_halfstarH[:,i]+F_halfstarH[:,i+1]-s*(Q_star[:,i+1]-Q_star[:,i])) # calculate F_halfstar with lower order Rusanov method
        
        F_halfstarL[:,0] = F[:,0] # make sure boundary conditions are applied
        F_halfstarL[:,ix-1] = F[:,ix-1]
        
        A_halfstar = F_halfstarH-F_halfstarL
        
        for i in range(1,ix-1): # Q_starstar
            Q_starstar[:,i] = Q[:,i]-(dt/dx)*(F_halfstarL[:,i]-F_halfstarL[:,i-1])
            
        Q_starstar[:,0] = Q[:,0] # Apply boundary conditions for Q_starstar
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        for i in range(1,ix-2): # Calculate corrected antidiffusion flux
            A_halfc[0,i] = np.sign(A_halfstar[0,i])*max(0,min(np.abs(A_halfstar[0,i]),np.sign(A_halfstar[0,i])*(Q_starstar[0,i+2]-Q_starstar[0,i+1])*dx/dt,np.sign(A_halfstar[0,i])*(Q_starstar[0,i]-Q_starstar[0,i-1])*dx/dt))
            A_halfc[1,i] = np.sign(A_halfstar[1,i])*max(0,min(np.abs(A_halfstar[1,i]),np.sign(A_halfstar[1,i])*(Q_starstar[1,i+2]-Q_starstar[1,i+1])*dx/dt,np.sign(A_halfstar[1,i])*(Q_starstar[1,i]-Q_starstar[1,i-1])*dx/dt))
            A_halfc[2,i] = np.sign(A_halfstar[2,i])*max(0,min(np.abs(A_halfstar[2,i]),np.sign(A_halfstar[2,i])*(Q_starstar[2,i+2]-Q_starstar[2,i+1])*dx/dt,np.sign(A_halfstar[2,i])*(Q_starstar[2,i]-Q_starstar[2,i-1])*dx/dt))
        
        A_halfc[:,0] = A_half[:,0] # Need to apply boundary conditions on the corrected antidiffusion flux
        A_halfc[:,ix-2] = A_half[:,ix-2]
        A_halfc[:,ix-1] = A_half[:,ix-1]
        
        for i in range(1,ix-1): # flux corrected solution for corrector
            Q_starstar[:,i] =  Q_starstar[:,i]-(dt/dx)*(A_halfc[:,i+1]-A_halfc[:,i]) 
       
        Q_starstar[:,0] = Q[:,0] # Apply boundary conditions for Q_starstar again just to be safe while I debug this
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        ### END OF CORRECTOR STEP
        
        Q = 0.5*(Q_star+Q_starstar) # Average Qstar and Qstarstar for next timelevel value of Q
        
        rho = Q[0,:] # Calculate all variables from Q
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
        
        F[0,:] = rho*u # Calculate new flux vector 
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a) # Advance timestep
        t = t+dt
        
    return rho, u, p, e, x

