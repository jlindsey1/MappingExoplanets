#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 21:12:41 2017

@author: Jordan
"""
print "Start..."

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import random as random
import matplotlib.ticker as mticker
from scipy.interpolate import interp1d
from sys import stdout

# Choose to include the middle region of the occultation
include_middle=True
if include_middle==True: inputfile='generated_data1'
if include_middle==False: inputfile='generated_data_nomid'
chi2file=np.genfromtxt(str(inputfile)+'.txt', names=True, delimiter=';',dtype=None)

# Comment out which model not needed. In this case, the no-limb-darkening model is commented out.
# This is the NULL HYPOTHESIS MODEL. <<<<<<<<<
modeldata1=np.genfromtxt('uniformingress1.txt', names=True, delimiter=';',dtype=None) #Uniform model
modeldata2=np.genfromtxt('uniformegress1.txt', names=True, delimiter=';',dtype=None)
#modeldata1=np.genfromtxt('nolimbingress1.txt', names=True, delimiter=';',dtype=None) #No-limb model
#modeldata2=np.genfromtxt('nolimbegress1.txt', names=True, delimiter=';',dtype=None)
graphspecs=np.genfromtxt('graph_specs.txt', names=True, delimiter=';',dtype=None)

# Import and define variables from the Jupyter Notebook.
P_total,P_full,P,flux_star,t_occultation,Initial,Length,Nslices=graphspecs['P_total'],graphspecs['P_full'],graphspecs['P'],graphspecs['flux_star'],graphspecs['t_occultation'],graphspecs['Initial'],graphspecs['Length'],graphspecs['Nslices']
print P_total,P_full,P,flux_star,t_occultation,Initial,Length,Nslices
P_total_initial=P_total
P_full_initial=P_full
Initial_Initial=Initial

# Define burn-in and main run:
Nburn=int(1e4)
Nsamples=int(1e5)

# Define initial step size and priors:
mid_step=0.0005
tot_step=0.0005
full_step=0.0005
mid_min=0.49
mid_max=0.51
# Note Ttot and Tfull is half the true value, for simplicity on the model generation.
tot_min=0.005
tot_max=0.010
full_min=0.005
full_max=0.010
# This is the SD per point from the jupyter notebook.
sigma_value=35*1e-6

# Choose to save output figures.
savefigures=False

# Define function to re-generate the model found in the Jupyter notebook lightcurve section.
def generate_model(full,tot,mid,verbose):
    Initial=mid
    P_full=full
    P_total=tot
    if verbose==True: print "Details: ", Initial, P_full, P_total, Length
    plotrange=np.linspace(-P_total+Initial,-P_full+Initial, num=Nslices)
    plotrange2=np.linspace(P_full+Initial,P_total+Initial, num=Nslices)
    stepdifference=np.abs(plotrange[0]-plotrange[1])
    rangedifference=np.abs(plotrange2[0]-plotrange[-1])
    Nsteps_needed=int(round(rangedifference/stepdifference))
    plotrange3=np.linspace(plotrange[-1]+stepdifference,plotrange2[0]-stepdifference,num=Nsteps_needed)
    uniform_curve_x,uniform_curve_y=[],[]
    
    total_amount = np.sum(modeldata1['bin_values'])
    
    for i in range(Nslices):
        total_amount = total_amount - modeldata1['bin_values'][i]
        fractional_flux = (total_amount+flux_star)/(flux_star)
        uniform_curve_x.append(plotrange[i])
        uniform_curve_y.append(fractional_flux)
    
    if include_middle==True:
        for i in range(len(plotrange3)):
            uniform_curve_x.append(plotrange3[i])
            uniform_curve_y.append(1.) 
    
    total_amount = 0
        
    for i in range(Nslices):
        total_amount = total_amount + modeldata2['bin_values'][Nslices-i-1]
        fractional_flux = (total_amount+flux_star)/(flux_star)
        uniform_curve_x.append(plotrange2[i])
        uniform_curve_y.append(fractional_flux)
    
    maxvalue=np.max(uniform_curve_y)
    uniform_curve_x.append(1)
    uniform_curve_y.append(maxvalue)
    uniform_curve_x.insert(0,0)
    uniform_curve_y.insert(0,maxvalue)
    
    # Returns the two curve arrays of data
    return uniform_curve_x,uniform_curve_y

# Progress bar, as used in tje Jupyter Notebook.
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    stdout.write('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    stdout.flush()
    if iteration == total: 
        print()
    
# Define function to calculate the chi-squared.
def chi2_calculation(inputfile,sigma_value,full,tot,mid):
    # Rejects values outside the set priors:
    if (tot_min < tot < tot_max) and (mid_min < mid < mid_max) and (full_min < full < full_max) and (tot>full):
        interpolation_datax,interpolation_dataf=generate_model(full,tot,mid,verbose=False)
        interpolation = interp1d(interpolation_datax, interpolation_dataf)
        chi2_list=[]
        for i in range(len(chi2file)):
            difference=chi2file['flux_values'][i]-interpolation(float(chi2file['x_values'][i]))
            if sigma_value!=0: chi2_list.append(difference**2/sigma_value**2)
            else: chi2_list.append(difference**2)
        return np.exp(-0.5*np.sum(chi2_list)),np.sum(chi2_list),np.sum(chi2_list)/(len(chi2file)-2)
    else:
        # Returns a zero probability if rejected.
        return 0,-100,-100

# Plots the initial lightcurve data.
interpolation_datax,interpolation_dataf=generate_model(full=P_full,tot=P_total,mid=Initial,verbose=True)
plt.plot(interpolation_datax,interpolation_dataf,c='k')
plt.scatter(chi2file['x_values'],chi2file['flux_values'],c='r',s=12,lw=0)#,zorder=2)
if sigma_value!=0: plt.errorbar(chi2file['x_values'],chi2file['flux_values'],yerr=sigma_value,c='#696969',lw=1,ls='none')
plt.xlim(0.48,0.52)
plt.ylim(np.min(chi2file['flux_values']),np.max(chi2file['flux_values']))
plt.xlabel('Phase')
plt.ylabel('$F(t)/F$')
plt.ticklabel_format(useOffset=False)
plt.show()

# Sets the lightcyrve plot ranges
plotrange=np.linspace(-P_total+Initial,-P_full+Initial, num=Nslices)
plotrange2=np.linspace(P_full+Initial,P_total+Initial, num=Nslices)

# Defines values to create lightcurve.
stepdifference=np.abs(plotrange[0]-plotrange[1])
rangedifference=np.abs(plotrange2[0]-plotrange[-1])
Nsteps_needed=int(round(rangedifference/stepdifference))
plotrange3=np.linspace(plotrange[-1]+stepdifference,plotrange2[0]-stepdifference,num=Nsteps_needed)

# Generates plot, similar to Jupyter Notebook lightcurve.
fig, (ax1, ax2) = plt.subplots(1, 2, sharex='col', sharey='row')
ax1.set_xlim(-P_total+Initial,-P_full+Initial)
ax2.set_xlim(P_full+Initial,P_total+Initial)
ax1.set_ylim(0.99999, 1.0001)
ax2.set_ylim(0.99999, 1.0001)
ax1.ticklabel_format(axis='y',style='sci',useOffset=False)
myLocator = mticker.MultipleLocator(0.0003)
myLocator2 = mticker.MultipleLocator(0.0003)
ax1.xaxis.set_major_locator(myLocator)
ax2.xaxis.set_major_locator(myLocator2)
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax1.tick_params(labeltop='off')
ax2.yaxis.tick_right()
plt.subplots_adjust(wspace=0.10)
ax1.set_xlabel('Phase')
ax1.set_ylabel('$F(t)/F$')
ax1.xaxis.set_label_coords(1.05, -0.1)
d = .015 
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d,1+d), (-d,+d), **kwargs)
ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)

print "Starting burn-in..."

# Initialises variables and histogram lists
mid_old=Initial
tot_old=P_total
full_old=P_full
prob_new=0
accepted=float(0)
acceptance2=[]
hist_values1=[]
hist_values2=[]
hist_values3=[]
hist_values4=[]
hist_values5=[]
thinning=False
interpolation_datax,interpolation_dataf=generate_model(full_old,tot_old,mid_old,verbose=False)
ax1.scatter(interpolation_datax,interpolation_dataf)
ax2.scatter(interpolation_datax,interpolation_dataf)
for i in range(Nburn):
    # For each loop, increment parameters by a random amount between their step sizes.
    tot_test=tot_old+random.uniform(-1, 1)*tot_step
    full_test=full_old+random.uniform(-1, 1)*full_step
    mid_test=mid_old+random.uniform(-1, 1)*mid_step
    # Sample the PDF for the new and old parameters
    prob,chi2,redchi2 = chi2_calculation(inputfile,sigma_value,full_test,tot_test,mid_test)
    prob_old,chi2_old,redchi2_old = chi2_calculation(inputfile,sigma_value,full_old,tot_old,mid_old)
    if thinning==False:
        pass
    elif thinning==True and i%10==0:
        # If thinning is enabled, save only each 10th value.
        hist_values1.append(mid_old)
        hist_values2.append(tot_old)
        hist_values3.append(chi2_old)
    # Test if the new position has a higher prob than the old. If yes, accept it
    if prob > prob_old:
        prob_new=1
        tot_old=tot_test
        full_old=full_test
        mid_old=mid_test
        # Increment acceptance counter 
        accepted+=1
    elif prob == 0: pass
    else:
        # If not, determine probability ratios. Then generate a random number between 
        # 0 and 1. If this random number falls below the ratio, accept the new position.
        # If not, do not accept it. This allows the MCMC to exit local maxima. MH algorithm.
        prob_new=prob/prob_old
        acceptance_prob=random.random()
        if acceptance_prob <= prob_new:
            tot_old=tot_test
            full_old=full_test
            mid_old=mid_test
            accepted+=1
        else: pass
        
    # For every 50 samples, determine if acceptance ratio is above or below 0.25
    # If above, increase step size to lower acceptance rate and vice versa.
    if i%50 == 0:
        acceptance=accepted/50
        acceptance2.append(acceptance)
        randomno1=random.uniform(0, 1)
        randomno2=random.uniform(0, 1)
        randomno3=random.uniform(0, 1)
        if acceptance > 0.25:
            tot_step = tot_step+0.00005*randomno1
            full_step = full_step+0.00005*randomno2
            mid_step = mid_step+0.00005*randomno3
        elif acceptance < 0.25:
            tot_step = tot_step-0.00005*randomno1
            full_step = full_step-0.00005*randomno2
            mid_step = mid_step-0.00005*randomno3
        accepted=0.0
    
    printProgressBar(i + 1, Nburn , prefix = 'Progress:', suffix = 'Complete', length = 50)
       
# Print variables for sanity check
print "Acceptance rate: ", np.mean(acceptance2)
print "tot_step: ", tot_step
print "full_step: ", full_step
print "mid_step: ", mid_step

interpolation_datax,interpolation_dataf=generate_model(full_old,tot_old,mid_old,verbose=True)
ax1.scatter(interpolation_datax, interpolation_dataf,c='k',s=10)
ax2.scatter(interpolation_datax, interpolation_dataf,c='k',s=10)

print "Starting main sampling..."

thinning=False
prob_new=0
accepted=float(0)

for i in range(Nsamples):
    # Save samples per loop to histograms
    if thinning==False:
        hist_values1.append(mid_old)
        hist_values2.append(tot_old)
        hist_values3.append(chi2_old)
        hist_values4.append(redchi2_old)
        hist_values5.append(full_old)
    elif thinning==True and i%10==0:
        hist_values1.append(mid_old)
        hist_values2.append(tot_old)
        hist_values3.append(chi2_old)
        hist_values4.append(redchi2_old)
        hist_values5.append(full_old)
    # For each loop, increment parameters by a random amount between their step sizes
    tot_test=tot_old+random.uniform(-1, 1)*tot_step
    full_test=full_old+random.uniform(-1, 1)*full_step
    mid_test=mid_old+random.uniform(-1, 1)*mid_step
    # Sample the posterior PDF for the new and old parameters
    prob,chi2,redchi2 = chi2_calculation(inputfile,sigma_value,full_test,tot_test,mid_test)
    prob_old,chi2_old,redchi2_old = chi2_calculation(inputfile,sigma_value,full_old,tot_old,mid_old)
    # Test if the new position has a higher prob than the old. If yes, accept it
    if prob > prob_old:
        prob_new=1
        tot_old=tot_test
        full_old=full_test
        mid_old=mid_test
    elif prob == 0: pass
    else:
        # If not, determine probability ratios. Then generate a random number between 
        # 0 and 1. If this random number falls below the ratio, accept the new position.
        # If not, do not accept it. This allows the MCMC to exit local maxima. MH algorithm.
        prob_new=prob/prob_old
        acceptance_prob=random.random()
        if acceptance_prob <= prob_new:
            tot_old=tot_test
            full_old=full_test
            mid_old=mid_test
        else: pass
    printProgressBar(i + 1, Nsamples, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
# Find the mean and SD for each algorithm.
mean=np.mean(hist_values1)
standard_dev=np.std(hist_values1)
mean2=np.mean(hist_values2)
standard_dev2=np.std(hist_values2)
mean3=np.mean(hist_values5)
standard_dev3=np.std(hist_values5)
print "mean: ", mean, "SD: ", standard_dev
print "mean2: ", mean2, "SD2: ", standard_dev2
print "mean3: ", mean3, "SD2: ", standard_dev3

# Save the results:

with open('hist_values1.txt', 'w') as file_handler:
    for item in hist_values1:
        file_handler.write("{}\n".format(item))
        
with open('hist_values2.txt', 'w') as file_handler:
    for item in hist_values2:
        file_handler.write("{}\n".format(item))
        
with open('hist_values3.txt', 'w') as file_handler:
    for item in hist_values3:
        file_handler.write("{}\n".format(item))
        
with open('hist_values4.txt', 'w') as file_handler:
    for item in hist_values4:
        file_handler.write("{}\n".format(item))
        
with open('hist_values5.txt', 'w') as file_handler:
    for item in hist_values5:
        file_handler.write("{}\n".format(item))


# Plot the initial results with the end of the burn-in and the mean of the main run
interpolation_datax,interpolation_dataf=generate_model(mean3,mean2,mean,verbose=True)
ax1.scatter(interpolation_datax, interpolation_dataf,c='c',s=11)
ax1.scatter(chi2file['x_values'],chi2file['flux_values'],c='b',s=8,lw=0)#,zorder=2)
ax1.errorbar(chi2file['x_values'],chi2file['flux_values'],yerr=sigma_value,c='#696969',lw=1,ls='none')#,zorder=1)
ax2.scatter(interpolation_datax, interpolation_dataf,c='c',s=11)
ax2.scatter(chi2file['x_values'],chi2file['flux_values'],c='b',s=8,lw=0)#,zorder=2)
ax2.errorbar(chi2file['x_values'],chi2file['flux_values'],yerr=sigma_value,c='#696969',lw=1,ls='none')
if savefigures==True: plt.savefig('MCMClightcurve.pdf')
ax2.legend(fontsize='medium',loc=(-0.40,+0.65))
plt.show()

interpolation_datax,interpolation_dataf=generate_model(mean3,mean2,mean,verbose=False)
plt.plot(interpolation_datax,interpolation_dataf)
plt.scatter(chi2file['x_values'],chi2file['flux_values'],c='b',s=8,lw=0)#,zorder=2)
if sigma_value!=0: plt.errorbar(chi2file['x_values'],chi2file['flux_values'],yerr=sigma_value,c='#696969',lw=1,ls='none')
plt.xlim(0.49,0.51)
plt.ylim(np.min(chi2file['flux_values']),np.max(chi2file['flux_values']))
plt.xlabel('Phase')
plt.ylabel('$F(t)/F$')
plt.show()

print "Done."

