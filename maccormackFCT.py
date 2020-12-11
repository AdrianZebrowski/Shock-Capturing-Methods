#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 13:39:21 2020

@author: adrianzebrowski
"""
def maccormackFCT(case,c_max,dx,gamma):
    import numpy as np 
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
    F_halfstar = np.zeros((3,ix))
    F_halfL = np.zeros((3,ix))
    F_halfH = np.zeros((3,ix))
    A_half = np.zeros((3,ix))
    A_halfc = np.zeros((3,ix))
    
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
        for i in range(0,ix-1): # F_half  
            s = np.maximum(np.abs(u[i])+a[i],np.abs(u[i+1])+a[i+1])
            F_halfL[:,i] = 0.5*(F[:,i]+F[:,i+1]-s*(Q[:,i+1]-Q[:,i])) # calculate F_halfL using Rusanov
            
        F_halfL[:,0] = F[:,0]
        F_halfL[:,ix-1] = F[:,ix-1]
            
        for i in range(0,ix-1):
            F_halfH[:,i] = F[:,i+1] # calculate F_halfH
            
        F_halfH[:,0] = F[:,0]
        F_halfH[:,ix-1] = F[:,ix-1]
        
        A_half = F_halfH - F_halfL
        
        for i in range(1,ix-1): # intermediate MP solution
            Q_star[:,i] = Q[:,i]-(dt/dx)*(F_halfL[:,i]-F_halfL[:,i-1])
            
        Q_star[:,0] = Q[:,0]
        Q_star[:,ix-1] = Q[:,ix-1]
        
        for i in range(1,ix-2): # corrected antidiffusion flux
            A_halfc[0,i] = np.sign(A_half[0,i])*max(0,min(np.abs(A_half[0,i]),np.sign(A_half[0,i])*(Q_star[0,i+2]-Q_star[0,i+1])*dx/dt,np.sign(A_half[0,i])*(Q_star[0,i]-Q_star[0,i-1])*dx/dt))
            A_halfc[1,i] = np.sign(A_half[1,i])*max(0,min(np.abs(A_half[1,i]),np.sign(A_half[1,i])*(Q_star[1,i+2]-Q_star[1,i+1])*dx/dt,np.sign(A_half[1,i])*(Q_star[1,i]-Q_star[1,i-1])*dx/dt))
            A_halfc[2,i] = np.sign(A_half[2,i])*max(0,min(np.abs(A_half[2,i]),np.sign(A_half[2,i])*(Q_star[2,i+2]-Q_star[2,i+1])*dx/dt,np.sign(A_half[2,i])*(Q_star[2,i]-Q_star[2,i-1])*dx/dt))
        
        for i in range(1,ix-1): # flux corrected solution
            Q_star[:,i] = Q_star[:,i]-(dt/dx)*(A_halfc[:,i]-A_halfc[:,i-1])
            
        #Q_star[:,0] = Q[:,0]
        #Q_star[:,ix-1] = Q[:,ix-1]
          
        rho_star = Q_star[0,:]
        u_star = Q_star[1,:]/rho_star
        e_star = Q_star[2,:]/rho_star-u_star**2/2
        p_star = e_star*(gamma-1)*rho_star
        a_star = np.sqrt(gamma*p_star/rho_star)
    #    
        F[0,:] = rho_star*u_star
        F[1,:] = rho_star*u_star**2+p_star
        F[2,:] = rho_star*u_star*(e_star+p_star/rho_star+u_star**2/2)
    #        
        for i in range(1,ix-1): # F_half  
            s_star = np.maximum(np.abs(u_star[i])+a_star[i],np.abs(u_star[i+1])+a_star[i+1])
            F_halfL[:,i] = 0.5*(F[:,i]+F[:,i+1]-s_star*(Q_star[:,i+1]-Q_star[:,i])) # calculate F_halfL
    #        
        F_halfL[:,0] = F[:,0]
        F_halfL[:,ix-1] = F[:,ix-1]
    #        
        for i in range(0,ix):
            F_halfH[:,i] = F[:,i] # calculate F_halfH
            
        F_halfH[:,0] = F[:,0]
        F_halfH[:,ix-1] = F[:,ix-1]
        A_half = F_halfH - F_halfL
        
        for i in range(1,ix-1): # intermediate MP solution
            Q_starstar[:,i] = Q[:,i]-(dt/dx)*(F_halfL[:,i]-F_halfL[:,i-1])
            
        Q_starstar[:,0] = Q[:,0]
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        for i in range(1,ix-2): # corrected antidiffusion flux
            A_halfc[0,i] = np.sign(A_half[0,i])*max(0,min(np.abs(A_half[0,i]),np.sign(A_half[0,i])*(Q_starstar[0,i+2]-Q_starstar[0,i+1])*dx/dt,np.sign(A_half[0,i])*(Q_starstar[0,i]-Q_starstar[0,i-1])*dx/dt))
            A_halfc[1,i] = np.sign(A_half[1,i])*max(0,min(np.abs(A_half[1,i]),np.sign(A_half[1,i])*(Q_starstar[1,i+2]-Q_starstar[1,i+1])*dx/dt,np.sign(A_half[1,i])*(Q_starstar[1,i]-Q_starstar[1,i-1])*dx/dt))
            A_halfc[2,i] = np.sign(A_half[2,i])*max(0,min(np.abs(A_half[2,i]),np.sign(A_half[2,i])*(Q_starstar[2,i+2]-Q_starstar[2,i+1])*dx/dt,np.sign(A_half[2,i])*(Q_starstar[2,i]-Q_starstar[2,i-1])*dx/dt))
        
        for i in range(1,ix-1): # flux corrected solution
            Q_starstar[:,i] = Q_starstar[:,i]-(dt/dx)*(A_halfc[:,i]-A_halfc[:,i-1])
            
        Q_starstar[:,0] = Q[:,0]
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        Q = 0.5*(Q_star+Q_starstar)
            
        #Q[:,:] = Q_star
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