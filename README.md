# GRS_voids

# GRS-voids
Using Fermi-LAT 4LAC catalog of gamma ray sources and SDSS DR 9/11 catalogs of optically detected quasars in combination with void catalog from Sutter et. al. to calculate fraction of void along line of sight to sources, investigations of observation that gamma-ray detected sources are along lines of sight with a larger void fraction. 

This README document outlines the catalogs and code necessary to reproduce results in both Furniss 2025 (10.3847/2041-8213/adae9d) and Jackson 202x (in prep.). Processed catalogs are included, as well as the code to process the raw catalogs. The order is as follows.

CONTAINS:
Folder - "Raw Data" - contains both intitial files and output files from Part 1 below. [NOTE: NEED TO FIGURE OUT HOW TO GET SDSS_QSOS.XLSX IN HERE, CURRENTLY TOO BIG]

PART 1: OBTAINING DATA + DATA PARSING

1. Obtain void catalog created from SDSS DR7, DR9, DR11 (TO FIX: Include link to Sutter Public Cosmic Void Catalog - website currently down) included in the folder here 'SDSS Catalog'. Uses the 'sky_positions_all' from sdss_dr7, sdss_dr9, sdss_dr10 subfolders from the COMOVING data for bright1, bright2, dim1, dim2, lrgbright, lrgdim (dr7) cmassbright, cmassmid, cmassdim (dr9), and cmass1/2/3 and lowz1/2/3/4 (dr10).

2.  Obtain Fermi-LAT 4LAC catalog (https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/892/105#/browse), specifically highlat.dat and lowlat.dat. The catalog used in our study is a combined file of highlat.dat and lowlat.dat in 'FINALCorrectedRedshift.xlsx' [file name kept for historical purposes].

   
3.Obtain SDSS catalog of optically-detected quasars (https://cdsarc.cds.unistra.fr/viz-bin/cat/VII/289), specifically from (https://cdsarc.cds.unistra.fr/ftp/VII/289/DR16Q_Superset_v3.fits), also included in this repository as 'DR16Q_Superset_v3.fits'. 


4. The three input files of data are 'FINALCorrectedRedshift.xlsx' [Fermi-LAT gamma-ray detected sources], 'DR16Q_Superset_v3.fits' [SDSS optically-detectecd sources], and 'processedsutter_voids.xlsx'. Use dataParser.ipynb to go from folder 'SDSS Catalog' to 'processedsutter_voids.xlsx' and to go from 'DR16Q_Superset_v3.fits' to 'sdss_qsos.xlsx' (will need to download 'DR16Q_Superset_v3.fits'). 'FINALCorrectedRedshifts.xlsx' is just a combined file of highlat.dat and lowlat.dat from 4LAC catalog and is provided.

PART 2: OBTAINING VOIDINESS OF SOURCES IN VOID FOOTPRINT

1. Create necessary conda environment in descartes_fix.py (making sure patch.py is in the same folder) to fix an issue with the descartes package.
2. run voidiness_calc.py [NOTE: this requires the additional scripts voidiness.py, perimeter.py, footprint_filter.py and custom_functions.py to be in the same folder]. This will save two .xlsx files to a new folder titled exported_dataFrames with the voidiness for the 4LAC catalog and the SDSS Qso catalog.

We now have the voidiness of all the sources. 

PART 3: COMPARING VOIDINESS OF POPULATIONS - RESULTS FROM FURNISS ET. AL. 2025

1. Running redshift_matched.py and Void_KS_z-matched.py to reproduce results from Furniss et. al. (2025). Requires the 4lac and SDSS qso files with the calculated voidiness. Redshift_matched.py creates and saves a fits file into exported_dataFrames with the redshift matched samples of SDSS Qso sources, which is a required input into Void_KS_z-matched.py.
2. NOTE: there is a bug with the anchor bin in redshift_matched.py . If you get output "creating redshift_matched population..." and it sits there without showing a graph for more than 15 seconds, you need to try a different anchor bin because it gets stuck in a while-loop.

   

