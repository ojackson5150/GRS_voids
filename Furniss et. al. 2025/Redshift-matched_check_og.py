#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 13:53:48 2024

@author: Olivier Hervet

description: creates a redshift-matched population of SDSS QSO sources saved
into a fits file.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from astropy.io import fits
cmap = plt.get_cmap("tab10")

params = {#'backend': 'ps',
      'axes.titlesize': 16,
      'axes.labelsize': 15,
      'xtick.labelsize': 15,
      'ytick.labelsize': 15,
      'legend.fontsize': 15}
plt.rcParams.update(params)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

print("read 4LAC...")
four_lac = pd.read_excel('exported_dataFrames/4lacsutter_w_voidiness_dup_drop_above_z0_1.xlsx')
z_4LAC = four_lac["z"]
#select only redshifts <= 0.4
z_4LAC_04 = []
for i in range(len(z_4LAC)):
    if (z_4LAC[i] >= 0.4):
        z_4LAC_04.append(z_4LAC[i])


print("read SDSS...")
SDSS = pd.read_excel('exported_dataFrames/sdsssutter_w_voidiness_dup_drop_above_z0_1.xlsx')
#class_SDSS = SDSS["AUTOCLASS_PQN"]
z_SDSS = SDSS["z"]
Voidiness_SDSS = SDSS["Voidiness"]

#Joseph method
# is_qsr = SDSS.AUTOCLASS_PQN == "QSO"
# near_mask = (SDSS.Z_QN >= 0.1) & (SDSS.Z_QN < 0.4)
# near_qso_mask = is_qsr & near_mask
# near_Voidiness = SDSS[near_qso_mask].Voidiness
#end of Joseph methof


z_SDSS_QSO_04 = []
Voidiness_SDSS_QSO_04 = []
for i in range(len(z_SDSS)):
    if 0.4 <= z_SDSS[i] < 0.7:
        z_SDSS_QSO_04.append(z_SDSS[i])
        Voidiness_SDSS_QSO_04.append(Voidiness_SDSS[i])

#print(min(Voidiness_SDSS_QSO_04))

bin_edges = np.linspace(0.4,0.7,13)

plt.figure("Redshift distribution")
hist_z_4LAC = np.histogram(z_4LAC_04,bins=bin_edges)
norm_hist_z_4LAC = hist_z_4LAC[0] / float(len(z_4LAC_04)) #histogram normalization
plt.stairs(norm_hist_z_4LAC,edges=bin_edges,lw=2,label="Fermi 4LAC")

hist_z_SDSS = np.histogram(z_SDSS_QSO_04,bins=bin_edges)
norm_hist_z_SDSS = hist_z_SDSS[0] / float(len(z_SDSS_QSO_04)) #histogram normalization
plt.stairs(norm_hist_z_SDSS,edges=bin_edges,lw=2,label="SDSS QSOs")
# plt.xlabel(r'Redshift')
# plt.ylabel(r'$\#$ of sources')
# plt.legend()


#the larger sample need to be carved out to match the smaller sample.
#a redshift bin needs to be selected as anchor in the larger sample
#typically the bin i with the smallest ratio of sources large_sample[i]/small_sample[i]
#for our case, this is the second redshift bin of the SDSS
anchor_4LAC = hist_z_4LAC[0][1]
bins_frac_4LAC = hist_z_4LAC[0]/anchor_4LAC

anchor_SDSS = hist_z_SDSS[0][1]
bins_frac_SDSS = hist_z_SDSS[0]/anchor_SDSS

#how many sources need to keep in SDSS for each bin?

nb_to_keep = np.round(anchor_SDSS*(bins_frac_4LAC)).astype(int)
#nb_to_keep = np.minimum(nb_to_keep, hist_z_SDSS[0])

def z_matching(z_SDSS_QSO_04, Voidiness_SDSS_QSO_04, norm_hist_z_SDSS, nb_to_keep, bin_edges):
    """
    Method to create a redshift-matched distribution of a large sample of sources based on
    the redshift distribution of a smaller sample

    Parameters
    ----------
    z_SDSS_QSO_04 : list
        list of redshifts SDSS.
    Voidiness_SDSS_QSO_04 : list
        list of voidiness SDSS.
    norm_hist_z_SDSS : list
        normalized binned list of SDSS redshifts.
    nb_to_keep : list
        sources need to keep in SDSS for each bin.
    bin_edges : list
        redshifts at bin edges.

    Returns
    -------
    z_matched_SDSS: list
    Voidiness_matched_SDSS: list

    """
    z_matched_SDSS = []
    Voidiness_matched_SDSS = []
    hist_z_matched_SDSS = [0] * norm_hist_z_SDSS
    #print(len(hist_z_matched_SDSS))
    z_SDSS_QSO_tmp = z_SDSS_QSO_04[:]
    Voidiness_SDSS_QSO_tmp = Voidiness_SDSS_QSO_04[:]
    #pick a random source with redshift (kind of acceptance/reject method)
    while sum(hist_z_matched_SDSS) < sum(nb_to_keep):
        index = np.random.randint(len(z_SDSS_QSO_tmp))
        for i in range(1,len(bin_edges)):
            if bin_edges[i-1] <= z_SDSS_QSO_tmp[index] < bin_edges[i] and hist_z_matched_SDSS[i-1] < nb_to_keep[i-1]:
                z_matched_SDSS.append(z_SDSS_QSO_tmp[index])
                Voidiness_matched_SDSS.append(Voidiness_SDSS_QSO_tmp[index])
                hist_z_matched_SDSS[i-1] += 1
                #remove indexes to not use them again
                del z_SDSS_QSO_tmp[index]
                del Voidiness_SDSS_QSO_tmp[index]
                break
                #print(hist_z_matched_SDSS)
    return(z_matched_SDSS, Voidiness_matched_SDSS)


print("building z-matched sample...")

z_matched_SDSS, Voidiness_matched_SDSS = z_matching(z_SDSS_QSO_04, Voidiness_SDSS_QSO_04, norm_hist_z_SDSS, nb_to_keep, bin_edges)

hist_z_matched_SDSS = np.histogram(z_matched_SDSS,bins=bin_edges)
norm_hist_z_matched_SDSS = hist_z_matched_SDSS[0] / float(len(z_matched_SDSS)) #histogram normalization
plt.stairs(norm_hist_z_matched_SDSS,edges=bin_edges,lw=2,ls="--",color=cmap(3),label="SDSS QSOs z-matched")
plt.xlabel(r'Redshift')
plt.ylabel(r'Normalized $\#$ of sources')
plt.legend()
plt.show()



#write 500 sample of z_matched SDSS in fits file
fits_file = "exported_dataFrames/SDSSsutter_04-07_z-matched.fits"

bin_edges_voids = np.linspace(0.1,0.9,19)
nbins = 499
all_hist = [0] *nbins
COL = [0] *nbins*2
plt.figure("Voidiness")
plt.title(r"$0.1 \leq z < 0.4$")
for i in range(nbins):
    z_matched_SDSS, Voidiness_matched_SDSS = z_matching(z_SDSS_QSO_04, Voidiness_SDSS_QSO_04, norm_hist_z_SDSS, nb_to_keep, bin_edges)
    COL[i*2] = fits.Column(name='Redshift'+str(i+1), unit= 'z',  format='D', array=z_matched_SDSS)
    COL[i*2+1] = fits.Column(name='Voidiness'+str(i+1),  format='D', array=Voidiness_matched_SDSS)
    aa = plt.hist(Voidiness_matched_SDSS,bins=bin_edges_voids,fill=False, lw=2, alpha = 0.1, color="0.6",histtype='step',zorder=0)
    all_hist[i] = aa[0]
median_voids = np.median(all_hist,axis=0)
plt.stairs(median_voids,edges=bin_edges_voids,lw=2,color=cmap(3),label="median z-matched SDSS",zorder=1)
plt.xlabel(r'Voidiness')
plt.ylabel(r'$\#$ of sources')
plt.legend()
plt.show()

#fill fits file
print("Writing fits file...")
hdu = fits.BinTableHDU.from_columns(COL)
hdu.writeto(fits_file, overwrite=True)
