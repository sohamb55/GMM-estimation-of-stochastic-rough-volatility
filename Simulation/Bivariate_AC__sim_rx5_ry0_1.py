# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 17:23:03 2023

@author: Laptop
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from MFDFA import MFDFA
from MFDFA import fgn
from numbers import Number
import math

from scipy.stats import norm

from scipy.optimize import fsolve

def generate_fOU(sim, time_step, delta_t, xi_true1, lambda_true1, nu_true1, H1,\
                xi_true2, lambda_true2, nu_true2, H2, rho_y):
       
    
    # Generate fractional Gaussian noise with H index
    #dB = (t_final ** H) * fgn(N = t_final*int(1/delta_t), H = H)
    dB1 = np.array([fgn(N = time_step.size, H = H1) for i in range(sim)])
    
    dB2 = rho_y * dB1 + np.sqrt(1 - rho_y **2) * np.array([fgn(N = time_step.size, H = H2) for i in range(sim)])
    
    var_yt1 = (nu_true1**2)/(2*lambda_true1**(2*H1))*(math.gamma(1+2*H1))
    var_yt2 = (nu_true2**2)/(2*lambda_true2**(2*H2))*(math.gamma(1+2*H2))
    
    eta1 = math.log(xi_true1) -0.5*var_yt1
    eta2 = math.log(xi_true2) -0.5*var_yt2

    # Initialise the array y
    y1,y2 = np.zeros([sim,time_step.size]),np.zeros([sim,time_step.size])
    
    # Give some small random initial conditions
    y1[:,0]=np.random.normal(loc = eta1, scale = var_yt1, size = 1) #/ 10
    y2[:,0]=np.random.normal(loc = eta2, scale = var_yt2, size = 1) #/ 10


    # Integrate the process
    for i in range(1, time_step.size):
        y1[:,i] = eta1 + (y1[:,i-1] - eta1) * math.exp(-lambda_true1* delta_t) + \
        nu_true1 * math.exp(-lambda_true1* delta_t/2) * dB1[:,i]
        
        y2[:,i] = eta2 + (y2[:,i-1] - eta2) * math.exp(-lambda_true2* delta_t) + \
        nu_true2 * math.exp(-lambda_true2* delta_t/2) * dB2[:,i]
        

    return y1,y2

def generate_Xt(y1,y2,rho_x,rho_y, sim, time_step, delta_t):
    
    # Initialise the array x    
    #x = np.zeros([sim,time_step.size])
    x1,x2 = np.zeros([sim,time_step.size]),np.zeros([sim,time_step.size])
    
    #Compute spot volatility

    spot_sigma1, spot_sigma2 = np.exp(0.5*y1), np.exp(0.5*y2)
    # Give some small random initial conditions
    dW1 = np.random.normal(0, np.sqrt(delta_t), (sim,time.size))
    dW2 = rho_x * dW1 + np.sqrt(1 - rho_x **2) * np.random.normal(0, np.sqrt(delta_t), (sim,time_step.size))
    # Integrate the process
        
    for i in range(1, time.size):
        x1[:,i] = x1[:,i-1] +  spot_sigma1[:,i-1] * dW1[:,i] + \
                                        np.sqrt(rho_y*spot_sigma1[:,i-1] *spot_sigma2[:,i-1])*dW2[:,i]
        
        x2[:,i] = x2[:,i-1] +  np.sqrt(rho_y*spot_sigma1[:,i-1] *spot_sigma2[:,i-1])*dW1[:,i] +\
                                            spot_sigma2[:,i-1] * dW2[:,i]
    return x1,x2

def generate_RV(x,sim,samp_freq):
    
    x_samp = x[:,::samp_freq]
    num_splits = int(time.size / N)
    #print(num_splits)

    x_samp_daily = np.array([np.split(x_samp[i], num_splits) for i in range(sim)])
    RV_sq = np.square(np.diff(x_samp_daily))
    RV_daily = np.array([np.sum(RV_sq[i],axis =1) for i in range(sim)])  
    
    #print(RV_daily.shape)
    return RV_daily

def generate_RCoV(x1,x2,sim,samp_freq):
    
    x1_samp = x1[:,::samp_freq]
    x2_samp = x2[:,::samp_freq]
    num_splits = int(time.size / N)
    #print(num_splits)

    x1_samp_daily = np.array([np.split(x1_samp[i], num_splits) for i in range(sim)])
    x2_samp_daily = np.array([np.split(x2_samp[i], num_splits) for i in range(sim)])
    
    RCoV_prod = np.square(np.diff(x1_samp_daily)*np.diff(x2_samp_daily))
    RCoV_daily = np.array([np.sum(RCoV_prod[i],axis =1) for i in range(sim)])  
    
    #print(RV_daily.shape)
    return RCoV_daily

def draw_spot_vol_paths(y1,y2,rho_x,rho_y,filename):
    
    fig = plt.figure(figsize=(12, 6))
    #plt.style.use('default')
    plt.plot(np.exp(0.5*y1),label ='$H ^{(1)}=0.05$',color='blue')
    plt.plot(np.exp(0.5*y2),label = '$H ^{(2)}=0.3$',color='goldenrod')
    #plt.axhline(y=0, label ='$E[Y_t ^{(1)}] = E[Y_t ^{(2)}] = 0$',color='orange', linestyle='dashed',linewidth=2.0)
    plt.xlabel('t',fontsize=30)
    plt.ylabel('$Y_t$',fontsize=30)
    plt.title(r'Sample paths of $\sigma_t ^{(1)},\sigma_t ^{(2)} $ with $\rho_x = %1.1f $, $\rho_y = %1.1f $'%(rho_x,rho_y),fontsize=30)
    plt.legend(fontsize=16,loc="upper right")
    plt.grid()
    
    plt.savefig(filename)
    plt.show()
    
def draw_asset_paths(x1,x2,rho_x,rho_y,filename):
    
    fig = plt.figure(figsize=(10, 4))
    #plt.style.use('default')
    plt.plot(x1,label ='$H ^{(1)}=0.05$',color='blue')
    plt.plot(x2,label = '$H ^{(2)}=0.3$',color='goldenrod')
    plt.axhline(y=0, label ='$E[X_t ^{(1)}] = E[X_t ^{(2)}] = 0$',color='orange', linestyle='dashed',linewidth=2.0)
    plt.xlabel('t',fontsize=30)
    plt.ylabel('$X_t$',fontsize=30)
    plt.title(r'Sample paths of $X_t ^{(1)},X_t ^{(2)}$ with $\rho_x = %1.1f $, $\rho_y = %1.1f $'%(rho_x,rho_y),fontsize=20)
    plt.legend(fontsize=16,loc="upper right")
    plt.grid()
    
    plt.savefig(filename)
    plt.show()
    
    

def draw_RV_paths(rv1,rv2,rho_x,rho_y,filename):
    
    fig = plt.figure(figsize=(12, 6))
    #plt.style.use('default')
    plt.plot(rv1,label ='$H ^{(1)}=0.05$',color='slateblue')
    plt.plot(rv2,label = '$H ^{(2)}=0.3$',color='goldenrod')
    #plt.axhline(y=0, label ='$E[X_t ^{(1)}] = E[X_t ^{(2)}] = 0$',color='orange', linestyle='dashed',linewidth=2.0)
    plt.xlabel('t',fontsize=30)
    plt.ylabel('$RV_t$',fontsize=30)
    plt.title(r'Sample paths of $RV_t ^{(1)},RV_t ^{(2)}$ with $\rho =%1.1f$, $\rho_y = %1.1f$'%(rho_x, rho_y),fontsize=20)
    plt.legend(fontsize=16,loc="upper right")
    plt.grid()
    
    plt.savefig(filename)
    plt.show()

def draw_RCoV_paths(rcov,rho_x,rho_y,filename):
    
    fig = plt.figure(figsize=(12, 6))
    #plt.style.use('default')
    plt.plot(rcov,color='royalblue')
    #plt.plot(rv2,label = '$RV_t ^{(2)}$',color='goldenrod')
    #plt.axhline(y=0, label ='$E[X_t ^{(1)}] = E[X_t ^{(2)}] = 0$',color='orange', linestyle='dashed',linewidth=2.0)
    plt.xlabel('t',fontsize=30)
    plt.ylabel('$RCoV_t$',fontsize=30)
    plt.title(r'Sample paths of $RCoV_t ^{(1)}, t \in [t_0,T]$ with $\rho =%1.1f$, $\rho_y = %1.1f$'%(rho_x, rho_y),fontsize=20)
    #plt.legend(fontsize=16,loc="upper right")
    plt.grid()
    
    plt.savefig(filename)
    plt.show()
# Initialization of simulation parameters

t_final = 1000 #250#4000#4000#500
N = 390 #23400

delta_t = 1/N # The desired timestep of integration

sim = 200 #Number of simulations

# time array of the process
time = np.linspace(0, t_final, t_final * int(1 / delta_t))
rho_x,rho_y = 0.5,0.0

Y1,Y2 = generate_fOU(sim=sim, time_step=time, delta_t=delta_t, xi_true1 = 0.0225,\
                         lambda_true1 = 0.005, nu_true1 = 1.25, H1 = 0.05,\
                           xi_true2 = 0.0225, lambda_true2 = 0.015, nu_true2 = 0.5,\
                         H2 = 0.3 ,rho_y = rho_y)


X1,X2 = generate_Xt(y1=Y1, y2=Y2, rho_x = rho_x,rho_y =rho_y, sim=sim, time_step=time,delta_t=delta_t)

RV1_daily = generate_RV(x=X1,sim=sim, samp_freq=5)
RV2_daily = generate_RV(x=X2,sim=sim, samp_freq=5)
RCoV_daily = generate_RCoV(x1=X1,x2=X2,sim=sim, samp_freq=5)

np.savetxt('RV1_rhox'+str(rho_x)+'_rhoy_'+str(rho_y)+'.csv', RV1_daily, delimiter=",")
np.savetxt('RV2_rhox'+str(rho_x)+'_rhoy_'+str(rho_y)+'.csv', RV2_daily, delimiter=",")
np.savetxt('RCoV_rhox'+str(rho_x)+'_rhoy_'+str(rho_y)+'.csv', RCoV_daily, delimiter=",")

draw_spot_vol_paths(y1= Y1[0],y2=Y2[0],rho_x = rho_x,rho_y =rho_y, filename='Y_path_rhox_'+str(rho_x)+\
                 '_rhoy_'+str(rho_y)+'.png')

draw_asset_paths(x1= X1[0],x2=X2[0],rho_x = rho_x,rho_y =rho_y, filename='X_path_rhox_'+str(rho_x)+\
                 '_rhoy_'+str(rho_y)+'.png')
draw_RV_paths(rv1= RV1_daily[0],rv2=RV2_daily[0],rho_x = rho_x,rho_y =rho_y, filename='RV_path_rhox_'+str(rho_x)+\
                 '_rhoy_'+str(rho_y)+'.png')
draw_RCoV_paths(rcov= RCoV_daily[0],rho_x = rho_x,rho_y =rho_y, filename='RCoV_path_rhox_'+str(rho_x)+\
                 '_rhoy_'+str(rho_y)+'.png')
    

asset_corr = [np.corrcoef(X1[i],X2[i])[0][1] for i in range(sim)]
vol_corr = [np.corrcoef(np.exp(0.5*Y1[i]),np.exp(0.5*Y2[i]))[0][1] for i in range(sim)]
    
# =============================================================================
# a = 0.5
# R = np.var(Y2[0])
# 
# xi=0.0225
# lam = 0.015
# nu = 0.5
# H=0.3
# func = lambda tau : R -  ((nu**2)/(2*(lam)**tau))*math.gamma(1+2*tau)#1+2*((1.0 - np.exp(-tau))/(1.0 - np.exp(-a*tau))) 
# 
# tau_initial_guess = 0.2
# tau_solution = fsolve(func, tau_initial_guess)
# print(tau_solution)
#  
# =============================================================================
from scipy.stats import gaussian_kde
data = np.random.normal(0,1,1000) # Generate Data
density1 = gaussian_kde((asset_corr-np.mean(asset_corr))/np.std(asset_corr))
density2 = gaussian_kde((vol_corr-np.mean(vol_corr))/np.std(vol_corr))
density3 = gaussian_kde(data)

x_vals = np.linspace(-4,4,1000) # Specifying the limits of our data
density1.covariance_factor = lambda : .5 #Smoothing parameter
density2.covariance_factor = lambda : .5 #Smoothing parameter
density3.covariance_factor = lambda : .5 #Smoothing parameter
 
density1._compute_covariance()
density2._compute_covariance()
density3._compute_covariance()
fig = plt.figure(figsize=(10, 8))
plt.plot(x_vals,density1(x_vals),label=' Asset correlation ',color='darkgoldenrod')
plt.plot(x_vals,density2(x_vals),linestyle="dashdot",label ='Volatility correlation',color='darkblue')
plt.plot(x_vals,density3(x_vals),label='N(0,1)',color='black')
#plt.plot(x_vals,gaussian_kde(np.random.normal(0,1,200)))
plt.title('Dependent volatility', fontsize=20)
plt.ylabel('density',fontsize=30)
plt.legend(fontsize=16,loc="upper right")
plt.grid()
plt.savefig('asset_cor_0.5_vol_cor_0.0.png')
plt.show()