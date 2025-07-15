
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
import healpy
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats

params = {#'backend': 'ps',
      'axes.labelsize': 18,
      'xtick.labelsize': 17,
      'ytick.labelsize': 18,
      'legend.fontsize': 16}
plt.rcParams.update(params)
plt.rc('text', usetex=True)


filename = "detthresh_P8R3_12years_PL22.fits"

# Open the FITS file
with fits.open(filename) as hdul:
    # 2D array in Galactic coordinates
    data = hdul[0].data
    header = hdul[0].header
    
wcs_Fermi = WCS(header)
    
plt.figure("Fermi sensitivity map", figsize=(10,5))
ax = plt.subplot(projection=wcs_Fermi)    
#ax = plt.subplot()  
plt.imshow(data, origin='lower', cmap='cividis', aspect='equal')
overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='0.5', ls='dotted')
plt.colorbar(location='bottom', shrink=0.7, pad=0.05, label=r"Fermi sensitivity [erg cm$^{-2}$ s$^{-1}$]")
#plt.tight_layout()

# #SDSS footprint boundaries
# RA_bound = [121.3, 121.3, 257.03, 257.03] * u.degree
# DEC_bound = [-7.49, 61.68, -7.49, 61.68] * u.degree

# DEC_lim = np.linspace(-7.49,61.68,100)
# RA_lim = np.linspace(121.3,257.03,500)
# aa =[]
# bb =[]
# for i in range(len(RA_lim)):
#     for j in range(len(DEC_lim)):
#         aa.append(RA_lim[i])
#         bb.append(DEC_lim[j])

# Bounds_celestial = SkyCoord(aa* u.degree, bb* u.degree, frame='icrs')
# # Transform to galactic coordinates
# Bounds_galactic = Bounds_celestial.transform_to('galactic')
# #convert to pixels
# pix_x, pix_y = wcs_Fermi.wcs_world2pix(Bounds_galactic.l, Bounds_galactic.b, 1)
# plt.scatter(pix_x, pix_y)


four_lac = pd.read_excel('exported_dataFrames/4lacsutter_w_voidiness_dup_drop_above_z0_1.xlsx')
#print(four_lac.columns)
FLac_celestial = SkyCoord(four_lac.RAdeg.values* u.degree, four_lac.DEdeg.values* u.degree, frame='icrs')
# Transform to galactic coordinates
FLac_galactic = FLac_celestial.transform_to('galactic')
#convert to pixels
pix_x, pix_y = wcs_Fermi.wcs_world2pix(FLac_galactic.l, FLac_galactic.b, 1)
plt.scatter(pix_x, pix_y,color ="g", label=r"4LAC selected sources")
plt.tight_layout()
plt.legend()
plt.show()

int_pix_x = np.round(pix_x).astype(int)
int_pix_y = np.round(pix_y).astype(int)

sensitivity_4lac =[]
Void = []
for i in range(len(pix_x)):
    sensitivity_4lac.append(data[int_pix_y[i],int_pix_x[i]])
    Void.append(four_lac.Voidiness[i])
        
plt.figure("Voidiness")
plt.plot(Void,sensitivity_4lac,"o")

def linear(x,a,b):
    return a*x+b

coeff,pcov = curve_fit(linear, Void, sensitivity_4lac)

plt.plot(Void, linear(np.array(Void),coeff[0],coeff[1]),color='0', label=f"f = ax+b\n a={coeff[0]:.2} $\pm$ {np.sqrt(pcov[0][0]):.2}\n b={coeff[1]:.2} $\pm$ {np.sqrt(pcov[1][1]):.2}")

#1sigma range
X = np.linspace(0,max(Void),200)
rng = np.random.RandomState(seed=30)
parameter_samples = rng.multivariate_normal(coeff, pcov, 100000)
 
realizations = np.array([linear(X, pars[0], pars[1]) for pars in parameter_samples])
q = 100 * stats.norm.cdf(-1)    #1 is the 1 sigma
y_low = np.percentile(realizations, q, axis=0)
q = 100 * stats.norm.cdf(1)     #1 is the 1 sigma
y_high = np.percentile(realizations, q, axis=0)
plt.fill_between(X, y_low, y_high,facecolor='0.5',edgecolor='0.5',alpha = 0.5)
plt.grid(True, which="both", axis = 'both', linestyle='--')
plt.xlabel(r"Voidiness")
plt.ylabel(r"Fermi sensitivity [erg cm-2 s-1]")
plt.legend()
plt.tight_layout()
plt.show()


#David want to know the flux distribution of 4LAC sources and compare with Fermi LAT sensitivity
#sensitivity map E = 100MeV to 100 GeV
#corresponding flux in 4FGL=Energy_Flux100
FourLAC_fluxes = four_lac.Energy_Flux100
#in log scale
FourLAC_fluxes = np.log10(FourLAC_fluxes)
plt.figure("Fuxes_vs_theshold")
plt.hist(FourLAC_fluxes,bins=30)
max_Fermi_threshold = np.log10(max(sensitivity_4lac))
min_Fermi_threshold = np.log10(min(sensitivity_4lac))
plt.axvline(max_Fermi_threshold,color="0",label="maximum Flux theshold")
plt.axvline(min_Fermi_threshold,color="r",label="minimum Flux theshold")
plt.xlabel("Flux [erg cm-2 s-1]")
plt.ylabel("4LAC sources")
plt.legend()
plt.tight_layout()
plt.show()

FourLAC_PL_index = four_lac.PL_Index
Low_PL_index = []
High_PL_index = []
Low_void = []
High_void=[]
for i in range(len(FourLAC_fluxes)):
    if FourLAC_fluxes[i] < max_Fermi_threshold:
        Low_PL_index.append(FourLAC_PL_index[i])
        Low_void.append(Void[i])
    else:
        High_PL_index.append(FourLAC_PL_index[i])
        High_void.append(Void[i])
