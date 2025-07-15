import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import os
import sys
import re
from astropy.table import QTable
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo # Used to calculate comoving dist
import math
import descartes

print('*********** Beginning voidiness analysis ***********')

#assumes use of processedsutter_voids.xlsx
voids = pd.read_excel('processedsutter_voids.xlsx')

#1. ------------ Running perimeter.py --------------
print("Running perimeter.py: determines a 2D footprint of the cosmic voids")
import perimeter #expects processedsutter_voids.xlsx to be in the same folder


#4. ------------ data_processing.ipynb -------------
#filtering by redshift 
from footprint_filter import filter_by_footprint

def filter_by_redshift(voids_df, data_df):
    too_close = data_df.z < min(voids_df.z)
    data_df = data_df[~too_close]

    too_far = data_df.z > max(voids_df.z)
    return data_df[~too_far]

def add_cmvd(data_df):
    """Adds the comoving distance to data_df"""
    cmvd = cosmo.comoving_distance(data_df['z']) # Comoving distance to void center
    data_df['cmvd_Mpc'] = cmvd.value # add it to data table    
    return data_df

def filter_by_z_ra_dec(data_df, voids, cmvd_add = True, footprint_points_fn = "footprint_points_sutter.xlsx"):
    
    if type(data_df) == str:
        data_df = pd.read_excel(data_df)
    else:
        assert isinstance(data_df, pd.DataFrame)
        
    data_df = filter_by_redshift(voids,data_df)

    if cmvd_add:
        data_df = add_cmvd(data_df)
    data_df = filter_by_footprint(data_df, footprint_points_fn)
    return data_df

print('Filtering 4LAC sources within footprint: Using \'FINALCorrectedRedshifts.xlsx\'.')

fourlac = filter_by_z_ra_dec('Initial Data/FINALCorrectedRedshifts.xlsx', voids)
fourlac.to_excel("exported_dataFrames/z_ra_dec_filtered_4lac_sutter.xlsx")
print("Filtered 4lac: saved filtered voids to \'exported_dataFrames/z_ra_dec_filtered_4lac_sutter.xlsx\'")

## Massive SDSS DR16 Catalog

sdss = pd.read_excel('sdss_qsos.xlsx')
print('Filtering SDSS sources within footprint: Using \'sdss_qsos.xlsx\'.')
sdss_dr16 = filter_by_z_ra_dec(sdss, voids)
sdss_dr16.to_excel("exported_dataFrames/z_ra_dec_filtered_sdss_sutter.xlsx")
print("Filtered sdss : saved filtered voids to \'exported_dataFrames/z_ra_dec_filtered_sdss_sutter.xlsx\'")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(test_df['RAdeg'],test_df['z'],test_df['DEdeg'], color = 'darkorange')
ax.scatter(voids['RAdeg'],voids['z'],voids['DEdeg'], color = 'cornflowerblue', label = "Sutter Voids", alpha = 0.1)
ax.scatter(sdss_dr16['RAdeg'], sdss_dr16['z'], sdss_dr16['DEdeg'], color = "gold", label = "SDSS", alpha = 0.5)
ax.scatter(fourlac['RAdeg'], fourlac['z'], fourlac['DEdeg'], color = "maroon", label = "4lac")
ax.set_ylabel('redshift')
ax.set_xlabel('RA')
ax.set_zlabel('DEC')
ax.set_ylim(0, 1.5)

plt.legend()
plt.show()

print("4lac sources:", len(fourlac))
print("SDSS sources:", len(sdss_dr16))

#5. ------------ voidiness.ipynb ------------
#voidy calculation 
from voidiness import voidy_analysis
print("Running void calculation on dataset: 4lac")

four_lac = voidy_analysis(voids,'exported_dataFrames/z_ra_dec_filtered_4lac_sutter.xlsx' )
#four_lac.to_excel('exported_dataFrames/4lacsutter_w_voidiness.xlsx', index=False)

print("Running void calculation on dataset: SDSS")

sdss_void = voidy_analysis(voids,'exported_dataFrames/z_ra_dec_filtered_sdss_sutter.xlsx' )
#sdss.to_excel('exported_dataFrames/sdsssutter_w_voidiness.xlsx', index=False)

#dropping any duplicated sources and ensuring sources are above z = 0.1

sdss_void = sdss_void.drop_duplicates(subset=['RAdeg', 'DEdeg'], keep='first')
four_lac = four_lac.drop_duplicates(subset=['RAdeg', 'DEdeg'], keep='first')

four_lac[four_lac.z >=0.1].to_excel("exported_dataFrames/4lacdes_w_voidiness_dup_drop_above_z0_1.xlsx",index=False)
sdss_void[sdss_void.z >= 0.1].to_excel("exported_dataFrames/sdssdes_w_voidiness_dup_drop_above_z0_1.xlsx", index=False)

print("**********************************\n")
print("* Finished Voidiness Calculation *\n")
print("*        Exit with Success       *\n")
print("**********************************\n")
