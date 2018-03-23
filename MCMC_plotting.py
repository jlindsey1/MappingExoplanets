#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 21:12:41 2017

@author: Jordan
"""
# The following plots several figures from the MCMC. Not all are relevant, but the code has been kept in this single file
# for simplicity.

print "Start..."

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import SkewedGaussianModel
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Import MCMC results
with open('hist_values1.txt') as f: hist_values1 = f.read().splitlines()
with open('hist_values2.txt') as f: hist_values2 = f.read().splitlines()
with open('hist_values3.txt') as f: hist_values3 = f.read().splitlines()
with open('hist_values4.txt') as f: hist_values4 = f.read().splitlines()
with open('hist_values5.txt') as f: hist_values5 = f.read().splitlines()
hist_values1=[float(i) for i in hist_values1]
hist_values2=[float(i) for i in hist_values2]
hist_values3=[float(i) for i in hist_values3]
hist_values4=[float(i) for i in hist_values4]
hist_values5=[float(i) for i in hist_values5]

# Double Ttot and Tfull as only half values were used in the MCMC (to simplify maths)
hist_values2=np.array(hist_values2)*2
hist_values5=np.array(hist_values5)*2

include_middle=True
if include_middle==True: inputfile='generated_data1'
if include_middle==False: inputfile='generated_data_nomid'
chi2file=np.genfromtxt(str(inputfile)+'.txt', names=True, delimiter=';',dtype=None)
modeldata1=np.genfromtxt('uniformingress1.txt', names=True, delimiter=';',dtype=None) #Uniform model
modeldata2=np.genfromtxt('uniformegress1.txt', names=True, delimiter=';',dtype=None)
#modeldata1=np.genfromtxt('nolimbingress1.txt', names=True, delimiter=';',dtype=None) #No-limb model
#modeldata2=np.genfromtxt('nolimbgress1.txt', names=True, delimiter=';',dtype=None)

# Import graph specifications
graphspecs=np.genfromtxt('graph_specs.txt', names=True, delimiter=';',dtype=None)
P_total,P_full,P,flux_star,t_occultation,Initial,Length,Nslices=graphspecs['P_total'],graphspecs['P_full'],graphspecs['P'],graphspecs['flux_star'],graphspecs['t_occultation'],graphspecs['Initial'],graphspecs['Length'],graphspecs['Nslices']
print P_total,P_full,P,flux_star,t_occultation,Initial,Length,Nslices
P_total_initial=P_total*2
P_full_initial=P_full*2
Initial_initial=Initial

savefigures=False
sigma_value=35*1e-6 #SD per point
    
mean=np.mean(hist_values1)
median=np.median(hist_values1)
standard_dev=np.std(hist_values1)
mean2=np.mean(hist_values2)
median2=np.median(hist_values2)
standard_dev2=np.std(hist_values2)
mean3=np.mean(hist_values5)
median3=np.median(hist_values5)
standard_dev3=np.std(hist_values5)
print "mean: ", mean, "SD: ", standard_dev, "Median: ", median
print "mean2: ", mean2, "SD2: ", standard_dev2, "Median2: ", median2
print "mean3: ", mean3, "SD3: ", standard_dev3, "Median3: ", median3

# Defines the model generation function
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
    return uniform_curve_x,uniform_curve_y

interpolation_datax,interpolation_dataf=generate_model(0.00730,0.0080,0.50035,verbose=False)
plt.plot(interpolation_datax,interpolation_dataf)
plt.scatter(chi2file['x_values'],chi2file['flux_values'],c='b',s=8,lw=0)#,zorder=2)
if sigma_value!=0: plt.errorbar(chi2file['x_values'],chi2file['flux_values'],yerr=sigma_value,c='#696969',lw=1,ls='none')
plt.xlim(0.47,0.53)
plt.ylim(np.min(chi2file['flux_values']),np.max(chi2file['flux_values']))
plt.xlabel('Phase')
plt.ylabel('$F(t)/F$')
if savefigures==True: plt.savefig('final-mcmc-lightcurve1.pdf')
plt.show()

heatmap, xedges, yedges = np.histogram2d(hist_values1, hist_values3, bins=(100,100),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
contourplot=ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')

ax2.axis('off')
ax1.hist(hist_values1,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step")
ax4.hist(hist_values3,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step", orientation="horizontal")
ax3.axis('tight')
ax3.ticklabel_format(useOffset=False)
#myLocator = mticker.MultipleLocator(0.00)
#ax3.xaxis.set_major_locator(myLocator)
ax3.set_xlabel('Midpoint Phase Position')
ax3.set_ylabel('Chi-Squared Value')
ax1.set_ylabel('PDF')
ax4.set_xlabel('PDF')
ax3.set_xlim(np.min(hist_values1),np.max(hist_values1))
ax3.set_ylim(np.min(hist_values3)*0.95,np.max(hist_values3))
if savefigures==True: plt.savefig('chisquared-corner1.pdf')
plt.show()

plt.hist2d(hist_values1,hist_values3, bins=100)
plt.xlabel('Midpoint Phase Position')
plt.ylabel('Chi-Squared')
if savefigures==True: plt.savefig('chisquared-hist1.pdf')
plt.show()

plt.hist2d(hist_values2,hist_values3, bins=100)
plt.xlabel('Total Duration Phase')
plt.ylabel('Chi-Squared')
if savefigures==True: plt.savefig('chisquared-hist2.pdf')
plt.show()

plt.hist2d(hist_values1,hist_values3, bins=200)
plt.xlabel('Midpoint Phase Position')
plt.ylabel('Chi-Squared')
if savefigures==True: plt.savefig('chisquared-hist3.pdf')
plt.show()

plt.hist2d(hist_values2,hist_values3, bins=200)
plt.xlabel('Total Duration Phase')
plt.ylabel('Chi-Squared')
if savefigures==True: plt.savefig('chisquared-hist4.pdf')
plt.show()

plt.hist2d(hist_values5,hist_values3, bins=200)
plt.xlabel('Full Duration Phase')
plt.ylabel('Chi-Squared')
if savefigures==True: plt.savefig('chisquared-hist5.pdf')
plt.show()

heatmap, xedges, yedges = np.histogram2d(hist_values2, hist_values3, bins=(100,100),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
contourplot=ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')

ax2.axis('off')
ax1.hist(hist_values2,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step")
ax4.hist(hist_values3,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step", orientation="horizontal")
ax3.axis('tight')
ax3.ticklabel_format(useOffset=False)
#myLocator = mticker.MultipleLocator(0.00)
#ax3.xaxis.set_major_locator(myLocator)
ax3.set_xlabel('Total Duration Phase')
ax3.set_ylabel('Chi-Squared Value')
ax1.set_ylabel('Marginalised Chi-Squared PDF')
ax4.set_xlabel('Marginalised Chi-Squared PDF')
ax3.set_xlim(np.min(hist_values2),np.max(hist_values2))
ax3.set_ylim(np.min(hist_values3)*0.95,np.max(hist_values3))
if savefigures==True: plt.savefig('chisquared-corner2.pdf')
plt.show()

y,x,_=plt.hist(hist_values1,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.axvline(x=Initial_initial,c='k',lw=2,label='Origin')
plt.xlabel('Midpoint Phase Position')
plt.ylabel('Marginalised Chi-Squared PDF')
plt.ylim(0,y.max()*(1.05))
plt.vlines(x=(mean), ymin=0, ymax=y.max()*(1.05), color='g', label='Mean')
plt.vlines(x=(mean-standard_dev), ymin=0, ymax=y.max()*(1.05), color='r', label='$\sigma_-$')
plt.vlines(x=(mean-standard_dev*2), ymin=0, ymax=y.max()*(1.05), color='m', label='$2\sigma_-$')
plt.vlines(x=(mean+standard_dev), ymin=0, ymax=y.max()*(1.05), color='b', label='$\sigma_+$')
plt.vlines(x=(mean+standard_dev*2), ymin=0, ymax=y.max()*(1.05), color='c', label='$2\sigma_+$')
plt.legend()
if savefigures==True: plt.savefig('PDF1-modified.pdf')
plt.show()

n_hist, b_hist, patches_hist = plt.hist(hist_values1,bins=200,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.hist(hist_values1,bins=200,normed=1,facecolor="black",edgecolor='None',alpha=0.1,label="PDF")
plt.xlabel('Midpoint Phase Position')
plt.ylabel('Normalised PDF')
if savefigures == True: plt.savefig('plottemp.pdf')
bin_max = np.where(n_hist == n_hist.max())
print "Mode:", b_hist[bin_max][0]

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
        plt.axvline(Initial_initial,c='r')
        plt.axvline(maxvalx,c='k')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1
    #plt.axvline(xvals[i],c='b')

print testcondition
if savefigures == True: plt.savefig('asymmetric1.pdf')
plt.show()
###

y,x,_=plt.hist(hist_values2,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.axvline(x=P_total_initial,c='k',lw=2,label='Origin')
plt.xlabel('Total Duration Phase')
plt.ylabel('Marginalised Chi-Squared PDF')
plt.ylim(0,y.max()*(1.05))
plt.vlines(x=(mean2), ymin=0, ymax=y.max()*(1.05), color='g', label='Mean')
plt.vlines(x=(mean2-standard_dev2), ymin=0, ymax=y.max()*(1.05), color='r', label='$\sigma_-$')
plt.vlines(x=(mean2-standard_dev2*2), ymin=0, ymax=y.max()*(1.05), color='m', label='$2\sigma_-$')
plt.vlines(x=(mean2+standard_dev2), ymin=0, ymax=y.max()*(1.05), color='b', label='$\sigma_+$')
plt.vlines(x=(mean2+standard_dev2*2), ymin=0, ymax=y.max()*(1.05), color='c', label='$2\sigma_+$')
plt.legend()
if savefigures==True: plt.savefig('PDF2-modified.pdf')
plt.show()

n_hist, b_hist, patches_hist = plt.hist(hist_values2,bins=200,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.hist(hist_values2,bins=200,normed=1,facecolor="black",edgecolor='None',alpha=0.1,label="PDF")
plt.xlabel('Total Occultation Duration')
plt.ylabel('Normalised PDF')
if savefigures == True: plt.savefig('plottemp2.pdf')
bin_max = np.where(n_hist == n_hist.max())
print "Mode:", b_hist[bin_max][0]

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
        plt.axvline(maxvalx,c='k')
        plt.axvline(P_total_initial,c='r')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1

print testcondition
if savefigures == True: plt.savefig('asymmetric2.pdf')
plt.show()
###

y,x,_=plt.hist(hist_values5,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.axvline(x=P_full_initial,c='k',lw=2,label='Origin')
plt.xlabel('Full Duration Phase')
plt.ylabel('Marginalised Chi-Squared PDF')
plt.ylim(0,y.max()*(1.05))
plt.vlines(x=(mean3), ymin=0, ymax=y.max()*(1.05), color='g', label='Mean')
plt.vlines(x=(mean3-standard_dev3), ymin=0, ymax=y.max()*(1.05), color='r', label='$\sigma_-$')
plt.vlines(x=(mean3-standard_dev3*2), ymin=0, ymax=y.max()*(1.05), color='m', label='$2\sigma_-$')
plt.vlines(x=(mean3+standard_dev3), ymin=0, ymax=y.max()*(1.05), color='b', label='$\sigma_+$')
plt.vlines(x=(mean3+standard_dev3*2), ymin=0, ymax=y.max()*(1.05), color='c', label='$2\sigma_+$')
plt.legend()
if savefigures==True: plt.savefig('PDF3-modified.pdf')
plt.show()

n_hist, b_hist, patches_hist = plt.hist(hist_values5,bins=200,normed=1,edgecolor="black",facecolor="black",histtype="step",label="PDF")
plt.hist(hist_values5,bins=200,normed=1,facecolor="black",edgecolor='None',alpha=0.1,label="PDF")
plt.xlabel('Full Occultation Duration')
plt.ylabel('Normalised PDF')
if savefigures == True: plt.savefig('plottemp3.pdf')
bin_max = np.where(n_hist == n_hist.max())
print "Mode:", b_hist[bin_max][0]

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
        plt.axvline(maxvalx,c='k')
        plt.axvline(P_full_initial,c='r')
        plt.axvline(newx1[-1],c='#505050')
        plt.axvline(newx2[-1],c='#505050')
        print "Lower: ", np.abs(maxvalx-newx1[-1])
        print "Upper: ", np.abs(maxvalx-newx2[-1])
        break
    else: i+=1

print testcondition
if savefigures == True: plt.savefig('asymmetric3.pdf')
plt.show()
###

xpoints1=np.linspace(0,len(hist_values1),num=len(hist_values1))
xpoints2=np.linspace(0,len(hist_values2),num=len(hist_values2))
plt.scatter(xpoints1,hist_values1,c='r',s=3)
plt.xlabel('Number of Samples')
plt.ylabel('Midpoint Phase Position')
if savefigures==True: plt.savefig('parameter-variation1.pdf')
plt.show()
plt.scatter(xpoints2,hist_values2,c='b',s=3)
plt.xlabel('Number of Samples')
plt.ylabel('Total Duration Phase')
if savefigures==True: plt.savefig('parameter-variation2.pdf')
plt.show()
plt.scatter(xpoints2,hist_values5,c='b',s=3)
plt.xlabel('Number of Samples')
plt.ylabel('Full Duration Phase')
if savefigures==True: plt.savefig('parameter-variation3.pdf')
plt.show()

plt.scatter(xpoints2,hist_values4,c='m',s=3)
plt.xlabel('Number of Samples')
plt.ylabel('Reduced Chi Squared')
if savefigures==True: plt.savefig('parameter-variation3.pdf')
plt.show()

heatmap, xedges, yedges = np.histogram2d(hist_values1, hist_values2, bins=(100,100),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
contourplot=ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')
axins1 = inset_axes(ax3,
                    width="5%",  
                    height="92.5%",  
                    loc=1)

plt.colorbar(contourplot, cax=axins1, orientation="vertical")
ax2.axis('off')
ax1.hist(hist_values1,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step")
ax4.hist(hist_values2,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step", orientation="horizontal")
ax3.axis('tight')
ax3.ticklabel_format(useOffset=False)
myLocator = mticker.MultipleLocator(0.0003)
ax3.xaxis.set_major_locator(myLocator)
ax3.set_xlabel('Midpoint Position')
ax3.set_ylabel('Total Duration')
ax1.set_ylabel('Marginalised PDF')
ax4.set_xlabel('Marginalised PDF')
ax3.set_xlim(np.min(hist_values1),np.max(hist_values1))
ax3.set_ylim(np.min(hist_values2),np.max(hist_values2))
if savefigures==True: plt.savefig('corner-modified.pdf')
plt.show()

heatmap, xedges, yedges = np.histogram2d(hist_values1, hist_values5, bins=(100,100),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
contourplot=ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')
axins1 = inset_axes(ax3,
                    width="5%",  
                    height="92.5%",  
                    loc=1)

plt.colorbar(contourplot, cax=axins1, orientation="vertical")#, ticks=[1, 2, 3])
#plt.colorbar(contourplot,ax=ax3)
ax2.axis('off')
ax1.hist(hist_values1,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step")
ax4.hist(hist_values5,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step", orientation="horizontal")
ax3.axis('tight')
ax3.ticklabel_format(useOffset=False)
myLocator = mticker.MultipleLocator(0.0003)
ax3.xaxis.set_major_locator(myLocator)
ax3.set_xlabel('Midpoint Position')
ax3.set_ylabel('Full Duration')
ax1.set_ylabel('Marginalised PDF')
ax4.set_xlabel('Marginalised PDF')
ax3.set_xlim(np.min(hist_values1),np.max(hist_values1))
ax3.set_ylim(np.min(hist_values5),np.max(hist_values5))
if savefigures==True: plt.savefig('corner-modified2.pdf')
plt.show()


heatmap, xedges, yedges = np.histogram2d(hist_values2, hist_values5, bins=(100,100),normed=True)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
contourplot=ax3.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')
axins1 = inset_axes(ax3,
                    width="5%",  
                    height="92.5%",  
                    loc=1)

plt.colorbar(contourplot, cax=axins1, orientation="vertical")
ax2.axis('off')
ax1.hist(hist_values2,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step")
ax4.hist(hist_values5,bins=100,normed=1,edgecolor="black",facecolor="black",histtype="step", orientation="horizontal")
ax3.axis('tight')
ax3.ticklabel_format(useOffset=False)
#myLocator = mticker.MultipleLocator(0.00)
#ax3.xaxis.set_major_locator(myLocator)
ax3.set_xlabel('Total Duration')
ax3.set_ylabel('Full Duration')
ax1.set_ylabel('Marginalised PDF')
ax4.set_xlabel('Marginalised PDF')
ax3.set_xlim(np.min(hist_values2),np.max(hist_values2))
ax3.set_ylim(np.min(hist_values5),np.max(hist_values5))
if savefigures==True: plt.savefig('corner-modified3.pdf')
plt.show()

########################################


print "Done."

