#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:09:46 2018

@author: Jordan
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import SkewedGaussianModel

# This code achieves the same as for Ttot and Tfull, but instead for delta T midpoint parameter.

# Define number of samples.
Nsamples=int(1e5)

pi=np.pi

# Import midpoint results from the MCMC
with open('hist_values1.txt') as f: hist_values1 = f.read().splitlines()
hist_values1=[float(i) for i in hist_values1]
n_hist, b_hist, patches_hist = plt.hist(hist_values1,bins=200,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.xlabel('Midpoint Phase Position')
plt.ylabel('Normalised PDF')
bin_max = np.where(n_hist == n_hist.max())
print "Mode:", b_hist[bin_max][0]
plt.show()

# Define errors and positions for 0,90,180,270,360 degrees from the jupyter notebook for eccentricity alone.
w=[0, 90, 180, 270, 360]
midpoint_errors=[48.200245862989654, 28.057535610547685, 48.200245642213964, 28.057535610547685, 48.200245862989654]
midpoint_positions=[65786.525001015951, 65269.618991999996, 64752.712982984034, 65269.618991999996, 65786.525001015951]
combined_errors=[]

# n=? selects the value of omega in the above array to be used. IE: n=2 uses w=180deg
n=1

# Perform MC sample for modpoint position
for i in range(Nsamples):
    random_no=np.random.normal(0,1)
    error=midpoint_errors[n]*random_no
    combined_errors.append(midpoint_positions[n]+error)

# Define period
P=1.51087081*(3600*24)
hist_values1=np.array(hist_values1)
combined_errors=np.array(combined_errors)

# Calculate the combines eccentricity and volcano errors
combined_errors2=combined_errors+(hist_values1*P-P/2)

# Plot the two histograms together
n_hist, b_hist, patches_hist = plt.hist(combined_errors,bins=200,normed=1,histtype='step',edgecolor='k',label='No Offset ($\omega =$ {})'.format(w[n]))
bin_heights, bin_borders, _ = n_hist, b_hist, patches_hist
bin_center = bin_borders[:-1] + np.diff(bin_borders) / 2
xvals, yvals = bin_center, bin_heights
model = SkewedGaussianModel()
params = model.guess(yvals, x=xvals)
result = model.fit(yvals, params, x=xvals)
print result.fit_report()
plt.plot(xvals, result.best_fit,c='#808080')

n_hist, b_hist, patches_hist = plt.hist(combined_errors2,bins=200,normed=1,histtype='step',edgecolor='#FF66FF',label='Added Offset')
bin_heights, bin_borders, _ = n_hist, b_hist, patches_hist
bin_center = bin_borders[:-1] + np.diff(bin_borders) / 2
xvals, yvals = bin_center, bin_heights
model = SkewedGaussianModel()
params = model.guess(yvals, x=xvals)
result = model.fit(yvals, params, x=xvals)
print result.fit_report()
plt.plot(xvals, result.best_fit,c='m')

# Use gaussian fit with the combined histogram:

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
            if (position2>1e3*position1) and (position1!=0): raise Exception("Corresponding position for probability=({}) not correctly found. E1".format(position1))
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
        #plt.axvline(maxvalx,c='g')
        plt.axvline(maxvalx,c='k')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1
    #plt.axvline(xvals[i],c='b')

print testcondition
###

# Plot the results:

plt.xlabel('Midpoint Position (sec)')
plt.ylabel('Normalised PDF')
plt.legend(loc='upper right')
#plt.savefig('midpoint_mcmc2.pdf')
plt.show()

print "Done."