#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 13:39:21 2020

@author: adrianzebrowski
"""
def maccormackFV(case,c_max,dx,gamma):
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
        for i in range(1,ix-1): # Q_update
            Q_star[:,i] = Q[:,i]-(dt/dx)*(F[:,i+1]-F[:,i]) #Fstar right in loop
            
        Q_star[:,0] = Q[:,0]
        Q_star[:,ix-1] = Q[:,ix-1]
        
        rho_star = Q_star[0,:]
        u_star = Q_star[1,:]/rho_star
        e_star = Q_star[2,:]/rho_star-u_star**2/2
        p_star = e_star*(gamma-1)*rho_star
    
        F_halfstar[0,:] = rho_star*u_star
        F_halfstar[1,:] = rho_star*u_star**2+p_star
        F_halfstar[2,:] = rho_star*u_star*(e_star+p_star/rho_star+u_star**2/2)
        
        for i in range(1,ix-1): # Q_starstar
            Q_starstar[:,i] = Q[:,i]-(dt/dx)*(F_halfstar[:,i]-F_halfstar[:,i-1])
            
        Q_starstar[:,0] = Q[:,0]
        Q_starstar[:,ix-1] = Q[:,ix-1]
        
        Q = 0.5*(Q_star+Q_starstar)
        
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