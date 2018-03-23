#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:09:46 2018

@author: Jordan
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import SkewedGaussianModel

# Number of monte-carlo samples to generate for Tocc and Ttot
Nsamples=int(1e5)

# Choose to save output figures
savefigures=False

# Define parameters of the system:

b=0.126
b_p=0.092
b_m=0.078

P=1.51087081*(3600*24)
sigma_P=0.60*1e-6*3600*24
P_p=sigma_P

pi=np.pi

R_Sun=6.957e8
R_star=0.117*R_Sun
sigma_Rs=0.0036*R_Sun

R_Earth=6.3781e6
R_P=1.086*R_Earth
sigma_Rp=0.035*R_P

scale_param=20.50
scale_param_p=0.16
scale_param_m=0.31

inc=89.65
inc_p=0.22
inc_m=0.27

e=0.00622
sigma_e=0.00058
#sigma_e=0.00304
e_p=sigma_e
w=0*pi/180
sigma_w=3.11*pi/180
#sigma_w=34.24*pi/180
w_p=sigma_w

#######

# Define arrays to hold random values

R_star_errors=[]
R_P_errors=[]
inc_errors=[]
b_errors=[]
P_errors=[]

e_errors=[]
w_errors=[]

scale_param_errors=[]

# Generate random sample for each parameter:

for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=sigma_Rs*random_no
    R_star_errors.append(R_star+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=sigma_Rp*random_no
    R_P_errors.append(R_P+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    if random_no>0: error=b_p*random_no
    else: error=b_m*random_no
    b_errors.append(b+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    if random_no>0: error=inc_p*random_no
    else: error=inc_m*random_no
    inc_errors.append(inc+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    if random_no>0: error=scale_param_p*random_no
    else: error=scale_param_m*random_no
    scale_param_errors.append(scale_param+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=P_p*random_no
    P_errors.append(P+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=w_p*random_no
    w_errors.append(w+error)
    
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=e_p*random_no
    e_errors.append(e+error)

# Convert each parameter to numpy array

P_errors=np.array(P_errors)
scale_param_errors=np.array(scale_param_errors)
b_errors=np.array(b_errors)
w_errors=np.array(w_errors)
e_errors=np.array(e_errors)
R_star_errors=np.array(R_star_errors)
R_P_errors=np.array(R_P_errors)
inc_errors=np.array(inc_errors)

# Determine final Ttot and Tfull array

k_errors=R_P_errors/R_star_errors
const=(1-e_errors**2)**0.5 / (1-e_errors*np.sin(w_errors))
b_errors_mod=scale_param_errors*np.cos(inc_errors*pi/180)*((1-e_errors**2)/(1-e_errors*np.sin(w_errors)))
T_total=(P_errors/pi) * np.arcsin(scale_param_errors**-1 * ((1+k_errors)**2-b_errors_mod**2)**0.5 / (np.sin(inc_errors*pi/180)) ) * const
T_full=(P_errors/pi) * np.arcsin(scale_param_errors**-1 * ((1-k_errors)**2-b_errors_mod**2)**0.5 / (np.sin(inc_errors*pi/180)) ) * const

# Sometimes remove extreme values in Ttot as can interfere with algorithm.
T_total2=[]
for i in range(len(T_total)):
    if 2025 < T_total[i] < 2350: T_total2.append(T_total[i])

plt.clf()
plt.hist(b_errors,bins=200,normed=1,histtype='step',edgecolor='k')
plt.hist(b_errors_mod,bins=200,normed=1,histtype='step',edgecolor='c')
plt.show()

n_hist, b_hist, patches_hist = plt.hist(T_total2,bins=200,normed=1,histtype='step',edgecolor='k')
plt.xlabel('Total Occultation Duration (sec)')
plt.ylabel('Normalised PDF')

### CONFIDENCE INTERVAL SELECTOR: ########################################
# Determine x,y values for histogram peaks
bin_heights, bin_borders, _ = n_hist, b_hist, patches_hist
bin_center = bin_borders[:-1] + np.diff(bin_borders) / 2
xvals, yvals = bin_center, bin_heights
model = SkewedGaussianModel()
params = model.guess(yvals, x=xvals)
result = model.fit(yvals, params, x=xvals)
print result.fit_report()
plt.plot(xvals, result.best_fit,c='c',lw=2)

#Mode Finder:
maxval=0
maxvalx=0
for i in range(len(xvals)):
    if result.best_fit[i]>maxval:
        maxval=result.best_fit[i]
        maxvalx=xvals[i]

print "Curve Mode:", maxvalx

area = np.trapz(result.best_fit, x=xvals)#, dx=5)
print "area =", area 


# Perform calculation of 68% confidence range
summation1=0
summation2=0
prev_highest=[0]
prev_highest_position=[1e9]
i=0
newx1=[]
newy1=[]
newx2=[]
newy2=[]
# This algorithm will (from the LHS) find positions of equal height on the left and right side of the curve
# until the central region contains 68% of the total area.
# This is achieved using the gaussian fit and the trapezium rule.
while i < len(xvals):
    position1=result.best_fit[i]
    newx1.append(xvals[i])
    newy1.append(position1)
    summation1=np.trapz(newy1,x=newx1)
    found = False
    for j in range(len(xvals)):
        loc=len(xvals)-1-j
        if loc==-1: raise Exception("Array error.")
        position2=result.best_fit[loc]
        if (position2>=position1) and (found==False) and (xvals[loc]<=prev_highest_position[-1]) and (position2 >= prev_highest[-1]):
            if (position2>1e5*position1) and (position1!=0): raise Exception("Corresponding position for probability=({}) not correctly found. E1".format(position1))
            found = True
            prev_highest.append(position2)
            prev_highest_position.append(xvals[loc])
        if j>=len(n_hist) and found==False:
            raise Exception("Corresponding position for probability=({}) not found. E2".format(position1))
        if found == True:
            newx2.append(xvals[loc])
            newy2.append(position2)
            break
    summation2=np.abs(np.trapz(newy2,x=newx2))
    testcondition=1-(summation1+summation2)
    if testcondition<0.69:
        plt.axvline(maxvalx,c='k')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1

print testcondition
if savefigures == True: plt.savefig('ttot1.pdf')
plt.show()
###

# Sometimes remove extreme values in Tfull as can interfere with algorithm.
T_full2=[]
for i in range(len(T_full)):
    if 1600 < T_full[i] < 2000: T_full2.append(T_full[i])

n_hist, b_hist, patches_hist = plt.hist(T_full2,bins=200,normed=1,histtype='step',edgecolor='k')
plt.xlabel('Full Occultation Duration (sec)')
plt.ylabel('Normalised PDF')

### CONFIDENCE INTERVAL SELECTOR: ########################################
bin_heights, bin_borders, _ = n_hist, b_hist, patches_hist
bin_center = bin_borders[:-1] + np.diff(bin_borders) / 2
xvals, yvals = bin_center, bin_heights
model = SkewedGaussianModel()
params = model.guess(yvals, x=xvals)
result = model.fit(yvals, params, x=xvals)
print result.fit_report()
plt.plot(xvals, result.best_fit,c='c',lw=2)

#Mode Finder:
maxval=0
maxvalx=0
for i in range(len(xvals)):
    if result.best_fit[i]>maxval:
        maxval=result.best_fit[i]
        maxvalx=xvals[i]

print "Curve Mode:", maxvalx

area = np.trapz(result.best_fit, x=xvals)#, dx=5)
print "area =", area 

summation1=0
summation2=0
prev_highest=[0]
prev_highest_position=[1e9]
i=0
newx1=[]
newy1=[]
newx2=[]
newy2=[]
# This algorithm will (from the LHS) find positions of equal height on the left and right side of the curve
# until the central region contains 68% of the total area.
# This is achieved using the gaussian fit and the trapezium rule.
while i < len(xvals):
    position1=result.best_fit[i]
    newx1.append(xvals[i])
    newy1.append(position1)
    summation1=np.trapz(newy1,x=newx1)
    found = False
    for j in range(len(xvals)):
        loc=len(xvals)-1-j
        if loc==-1: raise Exception("Array error.")
        position2=result.best_fit[loc]
        if (position2>=position1) and (found==False) and (xvals[loc]<=prev_highest_position[-1]) and (position2 >= prev_highest[-1]):
            if (position2>1e5*position1) and (position1!=0): raise Exception("Corresponding position for probability=({}) not correctly found. E1".format(position1))
            found = True
            prev_highest.append(position2)
            prev_highest_position.append(xvals[loc])
            #plt.axvline(xvals[loc],c='m')
        if j>=len(n_hist) and found==False:
            raise Exception("Corresponding position for probability=({}) not found. E2".format(position1))
        if found == True:
            newx2.append(xvals[loc])
            newy2.append(position2)
            break
    summation2=np.abs(np.trapz(newy2,x=newx2))
    testcondition=1-(summation1+summation2)
    if testcondition<0.69:
        plt.axvline(maxvalx,c='k')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1
    #plt.axvline(xvals[i],c='b')

print testcondition
if savefigures == True: plt.savefig('tfull1.pdf')
plt.show()
###

print "Done."