#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 10:30:27 2024

@author: Olivier Hervet
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pandas as pd
from scipy import stats
from scipy.special import erfinv
from astropy.io import fits
cmap = plt.get_cmap("tab10")

params = {#'backend': 'ps',
      'axes.labelsize': 15,
      'xtick.labelsize': 15,
      'ytick.labelsize': 15,
      'legend.fontsize': 15}
plt.rcParams.update(params)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def fix_hist_step_vertical_line_at_end(ax):
    axpolygons = [poly for poly in ax.get_children() if isinstance(poly, Polygon)]
    for poly in axpolygons:
        poly.set_xy(poly.get_xy()[:-1])


four_lac = pd.read_excel('exported_dataFrames/4lacsutter_w_voidiness_dup_drop_above_z0_1.xlsx')
z_4LAC = four_lac["z"]
#select only redshifts < 0.4
z_4LAC_04 = []
Voidiness = []
for i in range(len(z_4LAC)):
    if 0.4 <= z_4LAC[i] < 0.7:
        z_4LAC_04.append(z_4LAC[i])
        Voidiness.append(four_lac["Voidiness"][i])
Voidiness = np.array(Voidiness)

f = fits.open("exported_dataFrames/SDSSsutter_04-07_z-matched.fits")  # open a FITS file
tbdata = f[1].data  # assume the first extension is a table
tbcolumns = f[1].columns

nbins = 20
Pval_KS = [0]*499
stats_KS = [0]*499
Voidiness_SDSS = [0]*499
Z_matched_Voidiness_array = [0]*499
plt.figure("Voidiness z-matched samples")
plt.xlabel(r'Voidiness')
plt.ylabel(r'$\#$ of z-matched SDSS samples')
for i in range(499):
    plt.hist(tbdata["Voidiness"+str(i+1)],bins=nbins,fill=False,histtype='step',color='0.5',lw=2,alpha=0.1,label="Randomized coordinates")
    Voidiness_SDSS[i] = np.mean(tbdata["Voidiness"+str(i+1)])
    #KS test
    KS,pval = stats.ks_2samp(Voidiness, tbdata["Voidiness"+str(i+1)])
    #res = stats.anderson_ksamp([Voidiness, tbdata["Voidiness"+str(i+1)]], method=stats.PermutationMethod())
    Pval_KS[i] = pval
    stats_KS[i] = KS
    Z_matched_Voidiness_array[i] = tbdata["Voidiness"+str(i+1)]

Z_matched_Voidiness_array = np.array(Z_matched_Voidiness_array)


plt.figure("KS P-value distribution")
plt.hist(Pval_KS,bins=nbins,fill=False,histtype='step',color=cmap(0),lw=2)
plt.xlabel(r'P-value ')
plt.ylabel(r'$\#$ of z-matched SDSS + DES samples')
P_med = np.median(Pval_KS)
plt.axvline(P_med,color='r',label=f'median p-value: {P_med:.2e}')
plt.legend()

print("median P-value:", P_med)
print("Significance (one-tail):", np.sqrt(2) * erfinv(1-2*P_med), "Sigma")


plt.figure("KS stat distribution")
plt.hist(np.array(stats_KS),bins=nbins,fill=False,histtype='step',color=cmap(0),lw=2)
plt.xlabel(r'KS-value')
plt.ylabel(r'$\#$ of z-matched SDSS samples')
aa = np.median(stats_KS)
plt.axvline(aa,color='r',label=f'median KS-value: {aa:.4}')
plt.legend()
plt.show()



nbins = 2000
plt.figure("4LAC vs SDSS/DES z-matched CDF")
all_hist = [0]*499
for i in range(499):
    if i == 0:
        aa = plt.hist(tbdata["Voidiness"+str(i+1)],bins=nbins,cumulative=True,density=True, fill=False,histtype='step',color='0.5',lw=2,alpha=0.)
        bins = aa[1]
        all_hist[i] = list(aa[0])
    else:
        all_hist[i] = list(plt.hist(tbdata["Voidiness"+str(i+1)],bins=nbins,cumulative=True,density=True, fill=False,histtype='step',color='0.5',lw=2,alpha=0.)[0])


ax = plt.gca()
ax.hist(Voidiness,fill=False,bins=nbins,cumulative=True,density=True, histtype='step',color=cmap(3),label="Fermi 4LAC")


fix_hist_step_vertical_line_at_end(ax)

all_hist = np.array(all_hist)
median = np.zeros(nbins)
Onesigma_low = np.zeros(nbins)
Onesigma_up = np.zeros(nbins)
Twosigma_low = np.zeros(nbins)
Twosigma_up = np.zeros(nbins)
for i in range(nbins):
    median[i] = np.percentile(all_hist[:,i], 50)#np.median(all_hist[:,i])
    Onesigma_up[i] = np.percentile(all_hist[:,i], stats.norm.cdf(1)*100)
    Onesigma_low[i] = np.percentile(all_hist[:,i], stats.norm.cdf(-1)*100)
    Twosigma_up[i] = np.percentile(all_hist[:,i], stats.norm.cdf(2)*100)
    Twosigma_low[i] = np.percentile(all_hist[:,i], stats.norm.cdf(-2)*100)

#plt.hist(median,bins=nbins,density=True, fill=False,histtype='step',color='b')
binstep = (bins[0] + bins[1])
plt.plot(bins[:-1]+binstep/2, median, color="0", ds="steps-mid", label="median z-matched Sutter")
plt.plot(bins[:-1]+binstep/2, Onesigma_low, color="0.5", ds="steps-mid",label="1 sigma")
plt.plot(bins[:-1]+binstep/2, Onesigma_up, color="0.5", ds="steps-mid")
plt.plot(bins[:-1]+binstep/2, Twosigma_low, color="0.8", ds="steps-mid",label="2 sigma")
plt.plot(bins[:-1]+binstep/2, Twosigma_up, color="0.8", ds="steps-mid")
plt.xlabel(r'Voidiness')
plt.ylabel(r'CDF')
plt.legend(loc="upper left")
plt.title(r'$0.4 \leq z < 0.7$', fontsize = 14)
plt.show()

print("Mean SDSS Voidiness:", np.mean(Z_matched_Voidiness_array))

print("Mean Fermi Voidiness:", np.mean(Voidiness))




#remove 10 sources at random of the 4lac to check how it affects the p-values
# med_pval = [0]*100
# for j in range(100):
#     print(j)
#     z_test = z_4LAC_04[:]
#     voidiness_test = list(Voidiness[:])
#     for i in range(10):
#         index = np.random.randint(len(z_test))
#         del z_test[index]
#         del voidiness_test[index]
#     Pval_test = [0] *499
#     for k in range(499):
#         #KS test
#         KS,pval = stats.ks_2samp(voidiness_test, tbdata["Voidiness"+str(k+1)])
#         Pval_test[k] = pval
#     med_pval[j] = np.median(Pval_test)

# med_pval = np.array(med_pval)
# sigma_distrib = np.sqrt(2) * erfinv(1-2*med_pval)

# nbins = 15
# plt.figure("Sigma distrib")
# plt.hist(sigma_distrib,bins=nbins,fill=True,histtype='step')
# plt.xlabel("Significance KS test 4LAC vs SDSS zmatched")
# plt.ylabel("nb simulations with 10 sources removed 4LAC")

# std = np.std(sigma_distrib)

# plt.figure("KS P-value distribution")
# plt.hist(Pval_KS_sims,bins=nbins,fill=False,histtype='step',color=cmap(0),lw=2)
# plt.xlabel(r'P-value ')
# plt.ylabel(r'$\#$ of simulations')

# plt.figure("KS stat distribution * len(sample)")
# plt.hist(np.array(stats_KS_sims)*len(Voidiness),bins=nbins,fill=False,histtype='step',color=cmap(0),lw=2)
# plt.xlabel(r'KS-value *len(sample)')
# plt.ylabel(r'$\#$ of simulations')

# plt.figure("test")
# bin_edges = np.linspace(0,0.8,30)
# sim4 = plt.hist(four_lac["Sim_Voidiness_"+str(4)],cumulative=True,density=True,bins=bin_edges,fill=False,histtype='step',lw=2)
# sim6 = plt.hist(four_lac["Sim_Voidiness_"+str(6)],cumulative=True,density=True,bins=bin_edges,fill=False,histtype='step',lw=2)
# flac = plt.hist(Voidiness,fill=False,bins=bin_edges,cumulative=True,density=True, histtype='step',lw=2,color=cmap(3),label="Fermi 4LAC")
# print(stats.ks_2samp(Voidiness, four_lac["Sim_Voidiness_"+str(4)],method='exact'))
# print(stats.ks_2samp(Voidiness, four_lac["Sim_Voidiness_"+str(6)],method='exact'))

# sim4 = sim4[0]
# sim6 = sim6[0]
# flac = flac[0]
# sim4diff = 0
# sim6diff = 0
# #calculate the largest absolute difference in their dcf (should be sismiare to KS stats)
# for i in range(len(sim4)):
#     diff4 = np.abs(sim4[i]-flac[i])
#     if diff4 > sim4diff:
#         sim4diff = diff4
#     diff6 = np.abs(sim6[i]-flac[i])
#     if diff6> sim6diff:
#         sim6diff = diff6
