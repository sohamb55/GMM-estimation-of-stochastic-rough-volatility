# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 02:16:45 2023

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
import statsmodels.api as sm
from scipy.stats import norm

def generate_fOU(sim, time_step, delta_t, xi_true, lambda_true, nu_true, H):
    
    
    # Generate fractional Gaussian noise with H index
    #dB = (t_final ** H) * fgn(N = t_final*int(1/delta_t), H = H)
    dB = np.array([fgn(N = time_step.size, H = H) for i in range(sim)])
    var_yt = (nu_true**2)/(2*lambda_true**(2*H))*(math.gamma(1+2*H))
    eta = math.log(xi_true) -0.5*var_yt

    # Initialise the array y
    y = np.zeros([sim,time_step.size])
    #print(y.shape)
    # Give some small random initial conditions
    y[:,0]=np.random.normal(loc = eta, scale = var_yt, size = 1) #/ 10

    # Integrate the process
    for i in range(1, time_step.size):
        y[:,i] = eta + (y[:,i-1] - eta) * math.exp(-lambda_true* delta_t) + \
        nu_true * math.exp(-lambda_true* delta_t/2) * dB[:,i]

    return y

def generate_Xt(y,sim, time_step, delta_t):
    
    # Initialise the array x    
    x = np.zeros([sim,time_step.size])
    
    #Compute spot volatility
    spot_sigma = np.exp(0.5*y)
    # Give some small random initial conditions
    x[:,0]= 0 

    # Integrate the process
    for i in range(1, time_step.size):
        x[:,i] = x[:,i-1] +  spot_sigma[:,i-1] * np.sqrt(delta_t)*np.random.normal(loc = 0, scale = 1, size = 1)
    return x

def generate_RV(x,sim,samp_freq):
    
    x_samp = x[:,::samp_freq]
    num_splits = int(time.size / N)
    #print(num_splits)

    x_samp_daily = np.array([np.split(x_samp[i], num_splits) for i in range(sim)])
    RV_sq = np.square(np.diff(x_samp_daily))
    RV_daily = np.array([np.sum(RV_sq[i],axis =1) for i in range(sim)])  
    
    #print(RV_daily.shape)
    return RV_daily



def draw_fOU_paths(times, paths, sim, expectations, title=None, KDE=False, marginal=False, marginalT=None, envelope=False,
                lower=None, upper=None, style="seaborn-v0_8-whitegrid", colormap="RdYlBu_r", **fig_kw):
    with plt.style.context(style):
        if marginal:
            fig = plt.figure(**fig_kw)
            gs = GridSpec(1, 5)
            ax1 = fig.add_subplot(gs[:4])
            ax2 = fig.add_subplot(gs[4:], sharey=ax1)

            last_points = [path[-1] for path in paths]
            cm = plt.colormaps[colormap]
            n_bins = int(np.sqrt(sim))
            n, bins, patches = ax2.hist(last_points, n_bins, orientation='horizontal', density=True)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            col = bin_centers - min(bin_centers)  # scale values to interval [0,1]
            col /= max(col)
            for c, p in zip(col, patches):
                plt.setp(p, 'facecolor', cm(c))
            my_bins = pd.cut(last_points, bins=bins, labels=range(len(bins) - 1), include_lowest=True)
            colors = [col[b] for b in my_bins]

            if KDE:
                kde = sm.nonparametric.KDEUnivariate(last_points)
                kde.fit()  # Estimate the densities
                ax2.plot(kde.density, kde.support, '--', lw=1.75, alpha=0.6, label='$Y_t$  KDE', zorder=10)
                ax2.axhline(y=np.mean(last_points), linestyle='--', lw=1.75, label=r'$\overline{Y_t}$')
            else:
                marginaldist = marginalT
                x = np.linspace(marginaldist.ppf(0.005), marginaldist.ppf(0.995), 1000)
                ax2.plot(marginaldist.pdf(x), x, '-', lw=1.75, alpha=0.6, label='$Y_t$ pdf')
                ax2.axhline(y=marginaldist.mean(), linestyle='--', lw=1.75, label='$E[Y_t]$')

            plt.setp(ax2.get_yticklabels(), visible=False)
            ax2.set_title('$Y_t$')
            ax2.legend()

            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0, color=cm(colors[i]))
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[Y_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, alpha=0.25, color='grey')
            plt.subplots_adjust(wspace=0.025, hspace=0.025)

        else:
            fig, ax1 = plt.subplots(**fig_kw)
            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0)
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[Y_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, color='grey', alpha=0.25)

        fig.suptitle(title)
        #ax1.set_title('Simulated Paths $Y_t, t \in [t_0, T]$')  # Title
        ax1.set_xlabel('$t$',fontsize=20)
        ax1.set_ylabel('$Y_t$',fontsize=20)
        ax1.legend()
        #plt.show()

    return fig

def draw_asset_paths(times, paths, sim, expectations, title=None, KDE=False, marginal=False, marginalT=None, envelope=False,
                lower=None, upper=None, style="seaborn-v0_8-whitegrid", colormap="RdYlBu_r", **fig_kw):
    with plt.style.context(style):
        if marginal:
            fig = plt.figure(**fig_kw)
            gs = GridSpec(1, 5)
            ax1 = fig.add_subplot(gs[:4])
            ax2 = fig.add_subplot(gs[4:], sharey=ax1)

            last_points = [path[-1] for path in paths]
            cm = plt.colormaps[colormap]
            n_bins = int(np.sqrt(sim))
            n, bins, patches = ax2.hist(last_points, n_bins, orientation='horizontal', density=True)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            col = bin_centers - min(bin_centers)  # scale values to interval [0,1]
            col /= max(col)
            for c, p in zip(col, patches):
                plt.setp(p, 'facecolor', cm(c))
            my_bins = pd.cut(last_points, bins=bins, labels=range(len(bins) - 1), include_lowest=True)
            colors = [col[b] for b in my_bins]

            if KDE:
                kde = sm.nonparametric.KDEUnivariate(last_points)
                kde.fit()  # Estimate the densities
                ax2.plot(kde.density, kde.support, '--', lw=1.75, alpha=0.6, label='$X_t$  KDE', zorder=10)
                ax2.axhline(y=np.mean(last_points), linestyle='--', lw=1.75, label=r'$\overline{X_t}$')
            else:
                marginaldist = marginalT
                x = np.linspace(marginaldist.ppf(0.005), marginaldist.ppf(0.995), 1000)
                ax2.plot(marginaldist.pdf(x), x, '-', lw=1.75, alpha=0.6, label='$X_t$ pdf')
                ax2.axhline(y=marginaldist.mean(), linestyle='--', lw=1.75, label='$E[X_t]$')

            plt.setp(ax2.get_yticklabels(), visible=False)
            ax2.set_title('$Y_t$')
            ax2.legend()

            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0, color=cm(colors[i]))
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[X_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, alpha=0.25, color='grey')
            plt.subplots_adjust(wspace=0.025, hspace=0.025)

        else:
            fig, ax1 = plt.subplots(**fig_kw)
            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0)
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[X_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, color='grey', alpha=0.25)

        fig.suptitle(title)
        #ax1.set_title('Simulated Paths $Y_t, t \in [t_0, T]$')  # Title
        ax1.set_xlabel('$t$',fontsize=20)
        ax1.set_ylabel('$X_t$',fontsize=20)
        ax1.legend()
        #plt.show()

    return fig

# =============================================================================
# =============================================================================
# =============================================================================
def draw_RV_paths(times, paths, sim, expectations, title=None, KDE=False, marginal=False, marginalT=None, envelope=False,
               lower=None, upper=None, style="seaborn-v0_8-whitegrid", colormap="RdYlBu_r", **fig_kw):
    with plt.style.context(style):
        if marginal:
            fig = plt.figure(**fig_kw)
            gs = GridSpec(1, 5)
            ax1 = fig.add_subplot(gs[:4])
            ax2 = fig.add_subplot(gs[4:], sharey=ax1)

            last_points = [path[-1] for path in paths]
            cm = plt.colormaps[colormap]
            n_bins = int(np.sqrt(sim))
            n, bins, patches = ax2.hist(last_points, n_bins, orientation='horizontal', density=True)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            col = bin_centers - min(bin_centers)  # scale values to interval [0,1]
            col /= max(col)
            for c, p in zip(col, patches):
                plt.setp(p, 'facecolor', cm(c))
            my_bins = pd.cut(last_points, bins=bins, labels=range(len(bins) - 1), include_lowest=True)
            colors = [col[b] for b in my_bins]

            if KDE:
                kde = sm.nonparametric.KDEUnivariate(last_points)
                kde.fit()  # Estimate the densities
                ax2.plot(kde.density, kde.support, '--', lw=1.75, alpha=0.6, label='$RV_t$  KDE', zorder=10)
                ax2.axhline(y=np.mean(last_points), linestyle='--', lw=1.75, label=r'$\overline{RV_t}$')
            else:
                marginaldist = marginalT
                x = np.linspace(marginaldist.ppf(0.005), marginaldist.ppf(0.995), 1000)
                ax2.plot(marginaldist.pdf(x), x, '-', lw=1.75, alpha=0.6, label='$RV_t$ pdf')
                ax2.axhline(y=marginaldist.mean(), linestyle='--', lw=1.75, label='$E[RV_t]$')

            plt.setp(ax2.get_yticklabels(), visible=False)
            ax2.set_title('$Y_t$')
            ax2.legend()

            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0, color=cm(colors[i]))
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[RV_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, alpha=0.25, color='grey')
            plt.subplots_adjust(wspace=0.025, hspace=0.025)

        else:
            fig, ax1 = plt.subplots(**fig_kw)
            for i in range(sim):
                ax1.plot(times, paths[i], '-', lw=1.0)
            ax1.plot(times, expectations, '--', lw=1.75, label='$E[RV_t]$')
            if envelope:
                ax1.fill_between(times, upper, lower, color='grey', alpha=0.25)

        fig.suptitle(title)
        #ax1.set_title('Simulated Paths $Y_t, t \in [t_0, T]$')  # Title
        ax1.set_xlabel('$t$',fontsize=20)
        ax1.set_ylabel('$RV_t$',fontsize=20)
        ax1.legend()
        #plt.show()

    return fig


def univariate_sim(N,delta_t,sim,time,xi_true,lambda_true,nu_true,H_true):
    
    Y = generate_fOU(sim=sim, time_step=time, delta_t=delta_t, xi_true = xi_true,\
                       lambda_true = lambda_true, nu_true = nu_true, H = H_true)

    print('- Completion of fOU simulation')
    X = generate_Xt(y=Y, sim=sim, time_step=time,delta_t=delta_t)

    print('- Completion of log asset price simulation')

    RV_daily = generate_RV(x=X,sim=sim, samp_freq=5)
    #np.savetxt("RV7_daily_1.csv", RV_daily, delimiter=",")
    print('- Completion of Realized Variance simulation')


    expectations_Y7 = np.mean(Y,axis=0)
    marginal_dist_Y7 = norm(loc=np.mean(Y,axis=0)[-1], scale=np.std(Y[:,-1]))

# =============================================================================
#     fig_Y = draw_fOU_paths(times=time,paths=Y,sim=sim, expectations = expectations_Y7,\
#                      marginal=True,marginalT=marginal_dist_Y7,KDE=False,\
#                      title = 'Simulated paths $Y_t, t \in [0, T]$ for $H = {}$'.format(H_true),figsize=(12,6))
#     #fig_Y.savefig('fOU_h7_1.png')
# =============================================================================

    print('- Completion of fOU plot')

    expectations_X7 = np.mean(X,axis=0)
    marginal_dist_X7 = norm(loc=np.mean(X,axis=0)[-1], scale=np.std(X[:,-1]))

# =============================================================================
#     fig_X = draw_asset_paths(times=time,paths=X,sim=sim, expectations = expectations_X7,\
#                      marginal=True,marginalT=marginal_dist_X7,KDE=False,colormap='OrRd',\
#                      title = 'Simulated paths $X_t, t \in [0, T]$ for $H = {}$'.format(H_true),figsize=(12,6))
#     #fig.savefig('X_h7_1.png')
# =============================================================================

    print('- Completion of log asset price plot')

    expectations_RV7 = np.mean(RV_daily,axis=0)
    marginal_dist_RV7 = norm(loc=np.mean(RV_daily,axis=0)[-1], scale=np.std(RV_daily[:,-1]))

# =============================================================================
#     fig_RV = draw_RV_paths(times=np.arange(0,RV_daily.shape[1],1),paths=RV_daily,sim=sim, expectations = expectations_RV7,\
#                      marginal=True,marginalT=marginal_dist_RV7,KDE=False,\
#                      title = 'Simulated paths $RV_t, t \in [0, T]$ for $H = {}$'.format(H_true),figsize=(12,6))
#     #fig.savefig('RV_h7_1.png')
# =============================================================================

    print('- Completion of Realized Variance plot')

    return RV_daily#, fig_Y, fig_X, fig_RV



# Initialization of simulation parameters

t_final = 500 #250#4000#4000#500
N = 390 #23400

delta_t = 1/N # The desired timestep of integration

sim = 10 #Number of simulations

# time array of the process
time = np.linspace(0, t_final, t_final * int(1 / delta_t))

print("------------------------------------------------")
print('Begining of execution for H:',0.05)
print("------------------------------------------------")

RV05_daily= univariate_sim(N,delta_t,sim,time,xi_true=0.0225,\
                           lambda_true=0.005,nu_true=1.25,H_true=0.05)

np.savetxt("RV05_daily_1.csv", RV05_daily, delimiter=",")
fig_Y05.savefig('fOU_h05_1.png')
fig_X05.savefig('X_h05_1.png')
fig_RV05.savefig('RV_h05_1.png')

print("------------------------------------------------")
print('End of execution for H:',0.05)
print("------------------------------------------------")

###############################################################################

print("------------------------------------------------")
print('Begining of execution for H:',0.10)
print("------------------------------------------------")

RV1_daily = univariate_sim(N,delta_t,sim,time,xi_true=0.0225,\
                           lambda_true=0.01,nu_true=0.75,H_true=0.1)

np.savetxt("RV1_daily_1.csv", RV1_daily, delimiter=",")
fig_Y1.savefig('fOU_h1_1.png')
fig_X1.savefig('X_h1_1.png')
fig_RV1.savefig('RV_h1_1.png')

print("------------------------------------------------")
print('End of execution for H:',0.10)
print("------------------------------------------------")

###############################################################################


print("------------------------------------------------")
print('Begining of execution for H:',0.30)
print("------------------------------------------------")

RV3_daily = univariate_sim(N,delta_t,sim,time,xi_true=0.0225,\
                           lambda_true=0.015,nu_true=0.50,H_true=0.3)

np.savetxt("RV3_daily_1.csv", RV3_daily, delimiter=",")
fig_Y3.savefig('fOU_h3_1.png')
fig_X3.savefig('X_h3_1.png')
fig_RV3.savefig('RV_h3_1.png')

print("------------------------------------------------")
print('End of execution for H:',0.30)
print("------------------------------------------------")

###############################################################################

print("------------------------------------------------")
print('Begining of execution for H:',0.50)
print("------------------------------------------------")

RV5_daily= univariate_sim(N,delta_t,sim,time,xi_true=0.0225,\
                          lambda_true=0.0350,nu_true=0.3,H_true=0.5)

np.savetxt("RV5_daily_1.csv", RV5_daily, delimiter=",")
fig_Y5.savefig('fOU_h5_1.png')
fig_X5.savefig('X_h5_1.png')
fig_RV5.savefig('RV_h5_1.png')

print("------------------------------------------------")
print('End of execution for H:',0.50)
print("------------------------------------------------")

###############################################################################

print("------------------------------------------------")
print('Begining of execution for H:',0.7)
print("------------------------------------------------")

RV_7_daily = univariate_sim(N,delta_t,sim,time,xi_true=0.0225,\
                            lambda_true=0.07,nu_true=0.2,H_true=0.7)

np.savetxt("RV7_daily_1.csv", RV_7_daily, delimiter=",")
fig_Y7.savefig('fOU_h7_1.png')
fig_X7.savefig('X_h7_1.png')
fig_RV7.savefig('RV_h7_1.png')

print("------------------------------------------------")
print('End of execution for H:',0.7)
print("------------------------------------------------")