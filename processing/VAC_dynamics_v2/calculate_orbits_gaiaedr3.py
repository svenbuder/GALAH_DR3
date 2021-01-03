#!/usr/bin/env python
# coding: utf-8

# # Actions and Orbit caluclation with MC sampling for GALAH DR3 after Gaia eDR3
# 
# ## Author: Sven Buder
# 
# ### History:
# 201204 SB Created  
# 

# # What information you need
# 
# ra, dec, pmra, pmdec from Gaia eDR3  
# 
# distance:  
# if you want to use parallax: parallax and parallax_uncertainty  
# if you want to use covariances: covariance entries from Gaia eDR23 
# 
# vlos:  
# if you want to use rv_galah: rv_galah, e_rv_galah  
# if you want to use rv_gaia: rv_gaia, e_rv_gaia

# In[ ]:


# Preamble for notebook 

# Compatibility with Python 3
from __future__ import (absolute_import, division, print_function)

try:
    get_ipython().run_line_magic('matplotlib', 'inline')
    get_ipython().run_line_magic('config', "InlineBackend.figure_format='retina'")
except:
    pass

# Start timer
import time
start = time.time()

# Basic packages
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import os
import sys
import random
import glob
import pickle
import collections
import pandas

# Packages to work with FITS and (IDL) SME.out files
import astropy.io.fits as pyfits
import astropy.table as table
import astropy.coordinates as coord
import astropy.units as u
import math
from astropy.table import Table, hstack, vstack
from scipy.io.idl import readsav

# Matplotlib and associated packages for plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox,TransformedBbox
from matplotlib.image import BboxImage
from matplotlib.legend_handler import HandlerBase
from matplotlib._png import read_png
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors

params = {
    'font.family'        : 'sans',
    'font.size'          : 17,
    'axes.labelsize'     : 20,
    'ytick.labelsize'    : 16,
    'xtick.labelsize'    : 16,
    'legend.fontsize'    : 20,
    'text.usetex'        : True, 
    'text.latex.preamble': [r'\usepackage{upgreek}', r'\usepackage{amsmath}'],
    }   
plt.rcParams.update(params)

_parula_data = [[0.2081, 0.1663, 0.5292], 
                [0.2116238095, 0.1897809524, 0.5776761905], 
                [0.212252381, 0.2137714286, 0.6269714286], 
                [0.2081, 0.2386, 0.6770857143], 
                [0.1959047619, 0.2644571429, 0.7279], 
                [0.1707285714, 0.2919380952, 0.779247619], 
                [0.1252714286, 0.3242428571, 0.8302714286], 
                [0.0591333333, 0.3598333333, 0.8683333333], 
                [0.0116952381, 0.3875095238, 0.8819571429], 
                [0.0059571429, 0.4086142857, 0.8828428571], 
                [0.0165142857, 0.4266, 0.8786333333], 
                [0.032852381, 0.4430428571, 0.8719571429], 
                [0.0498142857, 0.4585714286, 0.8640571429], 
                [0.0629333333, 0.4736904762, 0.8554380952], 
                [0.0722666667, 0.4886666667, 0.8467], 
                [0.0779428571, 0.5039857143, 0.8383714286], 
                [0.079347619, 0.5200238095, 0.8311809524], 
                [0.0749428571, 0.5375428571, 0.8262714286], 
                [0.0640571429, 0.5569857143, 0.8239571429], 
                [0.0487714286, 0.5772238095, 0.8228285714], 
                [0.0343428571, 0.5965809524, 0.819852381], 
                [0.0265, 0.6137, 0.8135], 
                [0.0238904762, 0.6286619048, 0.8037619048], 
                [0.0230904762, 0.6417857143, 0.7912666667], 
                [0.0227714286, 0.6534857143, 0.7767571429], 
                [0.0266619048, 0.6641952381, 0.7607190476], 
                [0.0383714286, 0.6742714286, 0.743552381], 
                [0.0589714286, 0.6837571429, 0.7253857143], 
                [0.0843, 0.6928333333, 0.7061666667], 
                [0.1132952381, 0.7015, 0.6858571429], 
                [0.1452714286, 0.7097571429, 0.6646285714], 
                [0.1801333333, 0.7176571429, 0.6424333333], 
                [0.2178285714, 0.7250428571, 0.6192619048], 
                [0.2586428571, 0.7317142857, 0.5954285714], 
                [0.3021714286, 0.7376047619, 0.5711857143], 
                [0.3481666667, 0.7424333333, 0.5472666667], 
                [0.3952571429, 0.7459, 0.5244428571], 
                [0.4420095238, 0.7480809524, 0.5033142857], 
                [0.4871238095, 0.7490619048, 0.4839761905], 
                [0.5300285714, 0.7491142857, 0.4661142857], 
                [0.5708571429, 0.7485190476, 0.4493904762],
                [0.609852381, 0.7473142857, 0.4336857143], 
                [0.6473, 0.7456, 0.4188], 
                [0.6834190476, 0.7434761905, 0.4044333333], 
                [0.7184095238, 0.7411333333, 0.3904761905], 
                [0.7524857143, 0.7384, 0.3768142857], 
                [0.7858428571, 0.7355666667, 0.3632714286], 
                [0.8185047619, 0.7327333333, 0.3497904762], 
                [0.8506571429, 0.7299, 0.3360285714], 
                [0.8824333333, 0.7274333333, 0.3217], 
                [0.9139333333, 0.7257857143, 0.3062761905], 
                [0.9449571429, 0.7261142857, 0.2886428571], 
                [0.9738952381, 0.7313952381, 0.266647619], 
                [0.9937714286, 0.7454571429, 0.240347619], 
                [0.9990428571, 0.7653142857, 0.2164142857], 
                [0.9955333333, 0.7860571429, 0.196652381], 
                [0.988, 0.8066, 0.1793666667], 
                [0.9788571429, 0.8271428571, 0.1633142857], 
                [0.9697, 0.8481380952, 0.147452381], 
                [0.9625857143, 0.8705142857, 0.1309], 
                [0.9588714286, 0.8949, 0.1132428571], 
                [0.9598238095, 0.9218333333, 0.0948380952], 
                [0.9661, 0.9514428571, 0.0755333333], 
                [0.9763, 0.9831, 0.0538]]

parula = ListedColormap(_parula_data, name='parula')
parula_zero = _parula_data[0]
parula_0 = ListedColormap(_parula_data, name='parula_0')
parula_0.set_bad((1,1,1))
parula_r = ListedColormap(_parula_data[::-1], name='parula_r')

willi_blau = [0.0722666667, 0.4886666667, 0.8467]


# In[ ]:


debug = False


# ### Galpy initialization
# 
# We are using the McMillan17 potential from McMillan, 2017, MNRAS, 465, 76.  
# Contrary to galpy, its normalisation parameters are:  
# r_gc = 8.21 kpc (galpy: 8.0 kpc, Gravity Collaboration, 2018, A&A, 615, 15: 8.178 kpc).  
# v_gc = 233.1 km/s (galpy: 220 km/s)

# In[ ]:


import galpy
#from galpy.potential import MWPotential2014 as pot
from galpy.potential.mwpotentials import McMillan17 as pot
from galpy.util.bovy_conversion import get_physical
from galpy.actionAngle import actionAngleStaeckel
from galpy.orbit import Orbit

# Reference values
#r_galactic_centre = 8.178*u.kpc # Gravity Collaboration, 2019, A&A, 625, 10
r_galactic_centre = 8.21*u.kpc # McMillan Potential, 2017
z_galactic_plane = 25.0*u.pc # Bland-Hawthorn & Gerhard, 2016, ARA&A, 54, 529

print('Reference frame:')
print('R_GC = '+str(r_galactic_centre)+' (McMillan, 2017, MNRAS, 465, 76)')
print('phi_GC = '+str(0*u.rad))
print('z_GC = '+str(z_galactic_plane)+' (Bland-Hawthorn & Gerhard, 2016, ARA&A, 54, 529)')

v_total_sun = (np.tan(6.379*u.mas)*r_galactic_centre/u.yr).to(u.km/u.s) # pm_l by Reid & Brunthaler 2004, ApJ, 616, 872
print('V_total_sun: = '+"{:.2f}".format(v_total_sun)+' (Reid & Brunthaler 2004, ApJ, 616, 872)')
v_peculiar = [11.1, 15.17, 7.25]*u.km/u.s # U and W from Schoenrich, Binney, Dehnen, 2010, MNRAS, 403, 1829, V so that V = V_total-V_sun
print('V_peculiar = ',(v_peculiar),' (U and W from Schoenrich, Binney, Dehnen, 2010, MNRAS, 403, 1829)')
print('V-component of V_peculiar = 15.17 km/s, instead of 12.24 km/s by Schoenrich et al. (2010), for matching v_circular')
v_circular = np.round(v_total_sun-v_peculiar[1],1)
print('V_circular = ',(v_circular),' (McMillan, 2017, MNRAS, 465, 76)')

aAS = actionAngleStaeckel(
        pot   = pot,        #potential                                                                                                                                                                      
        delta = 0.45,       #focal length of confocal coordinate system                                                                                                                            
        c     = True        #use C code (for speed)                                                                                                                                                         
        )

#(RA = 17:45:37.224 h:m:s, Dec = −28:56:10.23 deg) (Reid& Brunthaler 2004)


# ### Let's get the Solar values

# In[ ]:


calculate_sun = True

if calculate_sun:
    sun = dict()

    # Create the Orbit instance
    o = Orbit(
        #ra, dec, dist, pm_ra, pm_dec, v_los
        vxvv=[0.*u.deg,0.*u.deg,0.*u.kpc,0.*u.mas/u.yr, 0.*u.mas/u.yr,0.*u.km/u.s],
        ro=r_galactic_centre,
        vo=v_circular,
        zo=z_galactic_plane,
        solarmotion=[-11.1, 15.17, 7.25]*u.km/u.s,
        #solarmotion='schoenrich',
        radec=True
    )

    #Galactocentric coordinates:
    sun['X_XYZ'] = o.helioX()#*u.kpc        
    sun['Y_XYZ'] = o.helioY()#*u.kpc
    sun['Z_XYZ'] = o.helioZ()#*u.kpc
    sun['U_UVW'] = o.U()#*u.km/u.s
    sun['V_UVW'] = o.V()#*u.km/u.s
    sun['W_UVW'] = o.W()#*u.km/u.s
    sun['R_Rzphi'] = o.R()#*u.kpc
    sun['phi_Rzphi'] = o.phi()#*u.rad
    sun['z_Rzphi'] = o.z()#*u.kpc
    sun['vR_Rzphi'] = o.vR()#*u.km/u.s
    sun['vT_Rzphi'] = o.vT()#*u.km/u.s        
    sun['vz_Rzphi'] = o.vz()#*u.km/u.s

    try:
        sun['J_R'], sun['L_Z'],sun['J_Z'], sun['omega_R'], sun['omega_phi'], sun['omega_z'], sun['angle_R'], sun['angle_phi'], sun['angle_z'] = aAS.actionsFreqsAngles(
            #R,vR,vT,z,vz[,phi]
            sun['R_Rzphi']*u.kpc,
            sun['vR_Rzphi']*u.km/u.s,
            sun['vT_Rzphi']*u.km/u.s,
            sun['z_Rzphi']*u.kpc,
            sun['vz_Rzphi']*u.km/u.s,
            sun['phi_Rzphi']*u.rad,
            ro=r_galactic_centre,vo=v_circular
        )
    except:
        sun['omega_R'] = [np.nan]
        sun['omega_phi'] = [np.nan]
        sun['omega_z'] = [np.nan]
        sun['angle_R'] = [np.nan]
        sun['angle_phi'] = [np.nan]
        sun['angle_z'] = [np.nan]
        try:
            sun['J_R'], sun['L_Z'],sun['J_Z'] = aAS(
            #R,vR,vT,z,vz[,phi]
            sun['R_Rzphi']*u.kpc,
            sun['vR_Rzphi']*u.km/u.s,
            sun['vT_Rzphi']*u.km/u.s,
            sun['z_Rzphi']*u.kpc,
            sun['vz_Rzphi']*u.km/u.s,
            sun['phi_Rzphi']*u.rad,
            ro=r_galactic_centre,vo=v_circular
        )
        except:
            sun['J_R'] = [np.nan]
            sun['L_Z'] = [np.nan]
            sun['J_Z'] = [np.nan]

    try:
        sun['ecc'], sun['zmax'], sun['R_peri'], sun['R_ap'] = aAS.EccZmaxRperiRap(
            #R,vR,vT,z,vz[,phi]
            sun['R_Rzphi']*u.kpc,
            sun['vR_Rzphi']*u.km/u.s,
            sun['vT_Rzphi']*u.km/u.s,
            sun['z_Rzphi']*u.kpc,
            sun['vz_Rzphi']*u.km/u.s,
            sun['phi_Rzphi']*u.rad,
            ro=r_galactic_centre,vo=v_circular
        )         
        sun['zmax']
        sun['R_peri']
        sun['R_peri']

    except:
        sun['ecc'] = [np.nan]
        sun['zmax'] = [np.nan]
        sun['R_peri'] = [np.nan]
        sun['R_ap'] = [np.nan]

    sun['Energy'] = o.E(pot=pot,ro=r_galactic_centre,vo=v_circular,zo=z_galactic_plane)

    print('Solar values:')
    print('X,Y,Z: '+"{:.2f}".format(sun['X_XYZ'])+' '+"{:.2f}".format(sun['Y_XYZ'])+' '+"{:.2f}".format(sun['Z_XYZ']))
    print('U,V,W: '+"{:.2f}".format(sun['U_UVW'])+' '+"{:.2f}".format(sun['V_UVW'])+' '+"{:.2f}".format(sun['W_UVW']))
    print('R,phi,z: '+"{:.2f}".format(sun['R_Rzphi'])+' '+"{:.2f}".format(sun['phi_Rzphi'])+' '+"{:.2f}".format(sun['z_Rzphi']))
    print('vR,vT,vz: '+"{:.2f}".format(sun['vR_Rzphi'])+' '+"{:.2f}".format(sun['vT_Rzphi'])+' '+"{:.2f}".format(sun['vz_Rzphi']))
    print('J_R,L_Z,J_Z: '+"{:.2f}".format(sun['J_R'][0])+' '+"{:.2f}".format(sun['L_Z'][0])+' '+"{:.2f}".format(sun['J_Z'][0]))
    print('Omega R/phi/z: '+"{:.2f}".format(sun['omega_R'][0])+' '+"{:.2f}".format(sun['omega_phi'][0])+' '+"{:.2f}".format(sun['omega_z'][0]))
    print('Angles R/phi/z: '+"{:.2f}".format(sun['angle_R'][0])+' '+"{:.2f}".format(sun['angle_phi'][0])+' '+"{:.2f}".format(sun['angle_z'][0]))
    print('ecc, zmax, R_peri, R_apo: '+"{:.2f}".format(sun['ecc'][0])+' '+"{:.2f}".format(sun['zmax'][0])+' '+"{:.2f}".format(sun['R_peri'][0])+' '+"{:.2f}".format(sun['R_ap'][0]))
    print('Energy: '+"{:.2f}".format(sun['Energy']))


# ### Input of 6D information in observable dimensions

# In[ ]:


try:
    galah_gaia_input = Table.read('/Users/svenbuder/GALAH_DR3/processing/VAC_dynamics_v2/VAC_Dynamics_V2_input.fits',1)
    out_dir = '/Users/svenbuder/GALAH_DR3/processing/VAC_dynamics_v2/'
except:
    galah_gaia_input = Table.read('/avatar/buder/trunk/GALAH/GALAH_DR3/processing/VAC_dynamics_v2/VAC_Dynamics_V2_input.fits',1)
    out_dir = '/avatar/buder/trunk/GALAH/GALAH_DR3/processing/VAC_dynamics_v2/'
        
full_length = len(galah_gaia_input['sobject_id'])
print("Initial nr. of entries")
print(full_length)

subset_size = 10000

try:
    subset = int(sys.argv[1])
except:
    subset = 0
    
plot_dynamics = False

if plot_dynamics:
    subset_size = 134
    subset = 13

if subset*subset_size >= full_length:
    sys.exit('The subset is beyond the length of GALAH DR3')

galah_gaia_input = galah_gaia_input[subset*subset_size:np.min([(subset+1)*subset_size,full_length])]

nr_galah_stars = len(galah_gaia_input['sobject_id'])
print("Nr. stars per subset")
print(nr_galah_stars)

nr_galah_stars_dynamics = np.where(
    (galah_gaia_input['use_dist_flag'] < 8) & # distance from BSTEP, BailerJones PhotoGeo/Geo, 1/Parallax
    (galah_gaia_input['use_rv_flag'] < 4) & # rv from Zwitter, GALAH, Gaia DR2
    np.isfinite(galah_gaia_input['pmra_pmdec_corr']) # Gaia eDR3 covariances available
)[0]

# This should only me activated for tests with subsets of GALAH DR3
#nr_galah_stars_dynamics = nr_galah_stars_dynamics[:100]

galah_gaia = galah_gaia_input[nr_galah_stars_dynamics]
nr_stars = len(galah_gaia['sobject_id'])
print("Selected number of stars")
print(nr_stars)


# In[ ]:


MC_size = 1000

print('MC Size: ',MC_size)

np.random.seed(123)

XYZ_labels       = ['X_XYZ','Y_XYZ','Z_XYZ']
UVW_labels       = ['U_UVW','V_UVW','W_UVW']

Rphiz_labels     = ['R_Rzphi','phi_Rzphi','z_Rzphi']
vRphiz_labels    = ['vR_Rzphi','vT_Rzphi','vz_Rzphi']

action_labels    = ['J_R','L_Z','J_Z']
frequency_labels    = ['omega_R','omega_phi','omega_z']
angle_labels    = ['angle_R','angle_phi','angle_z']
ext_orbit_labels = ['ecc', 'zmax', 'R_peri', 'R_ap', 'Energy']

orbit_labels = np.concatenate((
    XYZ_labels,
    UVW_labels,
    Rphiz_labels,
    vRphiz_labels,
    action_labels,
    frequency_labels,
    angle_labels,
    ext_orbit_labels    
    ))
print(orbit_labels)


# In[ ]:


# The final orbit information will go into a dictionary, which we initialise with np.nan values

orbit_information = collections.OrderedDict()

for each_orbit_label in orbit_labels:
    orbit_information[each_orbit_label] = np.zeros(nr_stars); orbit_information[each_orbit_label][:]=np.nan
for each_orbit_label in orbit_labels:
    orbit_information[each_orbit_label+'_5'] = np.zeros(nr_stars); orbit_information[each_orbit_label+'_5'][:]=np.nan
    orbit_information[each_orbit_label+'_50'] = np.zeros(nr_stars); orbit_information[each_orbit_label+'_50'][:]=np.nan
    orbit_information[each_orbit_label+'_95'] = np.zeros(nr_stars); orbit_information[each_orbit_label+'_95'][:]=np.nan


# In[ ]:


def distance2dmod(distance):
    """
    return the distance modulus for each distance
    """
    return np.array(5.0*np.log10(distance)-5.0)

def dmod2distance(dmod):
    """
    return the distance from distance modulus
    """
    return np.array(10**(1.0+0.2*dmod)).clip(min=0.0001)

def sample_distance_from_dmod(dist_16,dist_50,dist_84,MC_size):
    """
    Assume that the distance modulus is Gaussian
    (when using the mean the sigma from 50th-16th and 84th-16th percentile)
    """
    # Get distance modulus
    (dmod_16, dmod_50, dmod_84) = distance2dmod([dist_16,dist_50,dist_84])
    # Get mean sigma of distance modulus
    sigma_mean = 0.5*(dmod_84-dmod_16)

    # Sample with mean distance modulus and mean sigma of distance modulus
    dmod_sample = np.random.normal(dmod_50,scale=sigma_mean,size=MC_size)
    
    # return distance sample by converting distance modulus sample
    return(dmod2distance(dmod_sample))

def sample_distance_from_distance(dist_16,dist_50,dist_84,MC_size):
    """
    Assume that the distance modulus is Gaussian
    (when using the mean the sigma from 50th-16th and 84th-16th percentile)
    """
    # Get mean sigma of distance modulus
    sigma_mean = 0.5*(dist_84-dist_16)

    # return distance sample by sampling mean distance and mean sigma of distance
    return(np.random.normal(dist_50,scale=sigma_mean,size=MC_size).clip(min=0.0001))


def sample_distance_from2sidedGaussian(dist_16,dist_50,dist_84,MC_size):
    """
    Assume that the distance can be sampled via 2 Gaussians
    (one with sigma = 50th-16th and one with sigma 84th-50th)
    We sample them with MC_size around 0 and use their absolute values and then
    add/subtract from the 50th distance value.
    We then randomly select MC_size/2 of from either distribution.
    
    """
    distance_sigma_lo = np.abs(np.random.normal(loc = 0., scale = dist_50 -  dist_16, size = MC_size))
    distance_sigma_hi = np.abs(np.random.normal(loc = 0., scale = dist_84 -  dist_50, size = MC_size))
    select_lo_hi = np.concatenate((np.ones(int(MC_size/2)),np.zeros(MC_size-int(MC_size/2))))
    random.shuffle(select_lo_hi)
    
    return(np.array(
        (dist_50 + # 50th percentile
         select_lo_hi*distance_sigma_hi - # add selected sample from 84th-50th Gaussian
         (1-select_lo_hi)*distance_sigma_lo # subtract selected sample from 50th-16th Gaussian
        )).clip(min=0.0001) # Make sure distane is positive
    )

def sample_distance_from2sidedGaussian_DMod(dist_16,dist_50,dist_84,MC_size):
    """
    Assume that the distance can be sampled via 2 Gaussians
    (one with sigma = 50th-16th and one with sigma 84th-50th) in DMod
    We sample them with MC_size around 0 and use their absolute values and then
    add/subtract from the 50th distance value.
    We then randomly select MC_size/2 of from either distribution.
    
    """
    # Get distance modulus from distance
    (dmod_16, dmod_50, dmod_84) = distance2dmod([dist_16,dist_50,dist_84])

    # Sample Gaussian with sigmas according to uper and lower percentiles
    dmod_sigma_lo = np.abs(np.random.normal(loc = 0., scale = dmod_50 -  dmod_16, size = MC_size))
    dmod_sigma_hi = np.abs(np.random.normal(loc = 0., scale = dmod_84 -  dmod_50, size = MC_size))
    # Concatenate arrays with 1 (length = MC_size/2) and 0 (length = MC_size/2)
    select_lo_hi = np.concatenate((np.ones(int(MC_size/2)),np.zeros(MC_size-int(MC_size/2))))
    # Shuffle the values to randomly have randomly distributed values of 0 and 1
    random.shuffle(select_lo_hi)

    # Combine 2 different Gaussians sampled by the shuffled distribution
    dmod_sample = np.array(
        dmod_50 + # 50th percentile
        select_lo_hi*dmod_sigma_hi - # add selected sample from 84th-50th Gaussian
        (1-select_lo_hi)*dmod_sigma_lo # subtract selected sample from 50th-16th Gaussian
    )
    
    # return distance sample by converting distance modulus sample
    return((dmod2distance(dmod_sample)).clip(min=0.0001)) # Make sure distance is positive
    


# In[ ]:


def get_orbit_calculation_input(data, MC_size=1):
    """
    This function creates the 6D (ra, dec, distance, pmra, pmdec, vrad) input
    needed for galpy's Orbit() function and action calculations.    
    
    INPUT:
    data dictionary
    MC_size: if 1: only best value will be calculated
    
    OUTPUT:
    ra
    dec
    distance
    pmra
    pmdec
    vrad    
    """
    
    if MC_size==1:
        
        sample = dict()

        sample['ra'] = np.array([data['ra']])
        sample['dec'] = np.array([data['dec']])
        if data['use_dist_flag'] == 0:
            sample['distance'] = np.array([data['distance_bstep']])*1000.
        elif data['use_dist_flag'] == 1:
            sample['distance'] = np.array([data['r_med_photogeo']])
        elif data['use_dist_flag'] == 2:
            sample['distance'] = np.array([data['r_med_geo']])
        elif data['use_dist_flag'] == 4:
            sample['distance'] = np.array([1000./data['parallax_corr']])
        else:
            print("No useful distance available")
        sample['pmra'] = np.array([data['pmra']])
        sample['pmdec'] = np.array([data['pmdec']])
        if data['use_rv_flag'] == 0:
            sample['vrad'] = np.array([data['rv_obst']])
        elif data['use_rv_flag'] == 1:
            sample['vrad'] = np.array([data['rv_galah']])
        elif data['use_rv_flag'] == 2:
            sample['vrad'] = np.array([data['dr2_radial_velocity']])
        else:
            print("No useful radial velocity available")                                

    else:
        
        sample = dict()
        
        mu = np.array(
            [data['ra'],
             data['dec'],
             data['parallax_corr'],
             data['pmra'],
             data['pmdec'],
        ])
        
        s00 = (data['ra_error']/(3600.*1000.))**2
        s11 = (data['dec_error']/(3600.*1000.))**2
        s22 = data['parallax_error']**2
        s33 = data['pmra_error']**2
        s44 = data['pmdec_error']**2

        s01 = (data['ra_error']/(3600.*1000.)) * (data['dec_error']/(3600.*1000.)) * data['ra_dec_corr']
        s02 = (data['ra_error']/(3600.*1000.)) * data['parallax_error'] * data['ra_parallax_corr']
        s03 = (data['ra_error']/(3600.*1000.)) * data['pmra_error'] * data['ra_pmra_corr']
        s04 = (data['ra_error']/(3600.*1000.)) * data['pmdec_error'] * data['ra_pmdec_corr']

        s12 = (data['dec_error']/(3600.*1000.)) * data['parallax_error'] * data['dec_parallax_corr']
        s13 = (data['dec_error']/(3600.*1000.)) * data['pmra'] * data['dec_pmra_corr']
        s14 = (data['dec_error']/(3600.*1000.)) * data['pmdec'] * data['dec_pmdec_corr']

        s23 = data['parallax_error'] * data['pmra_error'] * data['parallax_pmra_corr']
        s24 = data['parallax_error'] * data['pmdec_error'] * data['parallax_pmdec_corr']

        s34 = data['pmra_error'] * data['pmdec_error'] * data['pmra_pmdec_corr']

        sigma = np.array([
            [s00, s01, s02, s03, s04],
            [s01, s11, s12, s13, s14],
            [s02, s12, s22, s23, s24],
            [s03, s13, s23, s33, s34],
            [s04, s14, s24, s34, s44]
        ])
        
        mu_sigma_sample = np.array(np.random.multivariate_normal(mu, sigma, size= MC_size))

        sample['ra'] = mu_sigma_sample[:,0]
        sample['dec'] = mu_sigma_sample[:,1]
        
        if data['use_dist_flag'] == 0:
            sample['distance'] = sample_distance_from2sidedGaussian_DMod(
                data['e16_distance_bstep']*1000,
                data['e50_distance_bstep']*1000,
                data['e84_distance_bstep']*1000,
                MC_size=MC_size
            )
        elif data['use_dist_flag'] == 1:
            sample['distance'] = sample_distance_from2sidedGaussian_DMod(
                data['r_lo_photogeo'],
                data['r_med_photogeo'],
                data['r_hi_photogeo'],
                MC_size=MC_size
            )
        elif data['use_dist_flag'] == 2:
            sample['distance'] = sample_distance_from2sidedGaussian_DMod(
                data['r_lo_geo'],
                data['r_med_geo'],
                data['r_hi_geo'],
                MC_size=MC_size
            )
        elif data['use_dist_flag'] == 4:
            sample['distance'] = (1000./mu_sigma_sample[:,2]).clip(min=0.0001)
        else:
            print("No useful distance available")

        sample['pmra'] = mu_sigma_sample[:,3]
        sample['pmdec'] = mu_sigma_sample[:,4]

        if data['use_rv_flag'] == 0:
            sample['vrad'] = np.array(np.random.normal(data['rv_obst'], scale=data['e_rv_obst'], size=MC_size))
        elif data['use_rv_flag'] == 1:
            sample['vrad'] = np.array(np.random.normal(data['rv_sme_v2'], scale=data['e_rv_sme'], size=MC_size))
        elif data['use_rv_flag'] == 2:
            sample['vrad'] = np.array(np.random.normal(data['dr2_radial_velocity'], scale=data['dr2_radial_velocity_error'], size=MC_size))
        else:
            print("No useful radial velocity available")
            
    return(sample)


# In[ ]:


def estimate_orbit_parameters(star_index, orbit_calculation_input):
    """
    Estimate orbit parameters from the given
    MC sample of 6D information for the Nr of stars
    and save it into orbit_information
    """

    # We are creating a dictionary for each star
    star_i = dict()

    ra     = orbit_calculation_input['ra']           *u.deg
    dec    = orbit_calculation_input['dec']          *u.deg
    dist   = orbit_calculation_input['distance']     *u.pc
    pm_ra  = orbit_calculation_input['pmra']         *u.mas/u.year
    pm_dec = orbit_calculation_input['pmdec']        *u.mas/u.year
    v_los  = orbit_calculation_input['vrad']         *u.km/u.s

    # Create the Orbit instance
    o = Orbit(
        vxvv=[ra,dec,dist,pm_ra, pm_dec,v_los],
        ro=r_galactic_centre,
        vo=v_circular,
        zo=z_galactic_plane,
        solarmotion=[-11.1, 15.17, 7.25]*u.km/u.s,
        radec=True
    )
    
    star_i = dict()
    #Galactocentric coordinates:
    star_i['X_XYZ'] = o.helioX()#*u.kpc        
    star_i['Y_XYZ'] = o.helioY()#*u.kpc
    star_i['Z_XYZ'] = o.helioZ()#*u.kpc
    star_i['U_UVW'] = o.U()#*u.km/u.s
    star_i['V_UVW'] = o.V()#*u.km/u.s
    star_i['W_UVW'] = o.W()#*u.km/u.s
    star_i['R_Rzphi'] = o.R()#*u.kpc
    star_i['phi_Rzphi'] = o.phi()#*u.rad
    star_i['z_Rzphi'] = o.z()#*u.kpc
    star_i['vR_Rzphi'] = o.vR()#*u.km/u.s
    star_i['vT_Rzphi'] = o.vT()#*u.km/u.s        
    star_i['vz_Rzphi'] = o.vz()#*u.km/u.s

    try:
        star_i['J_R'], star_i['L_Z'],star_i['J_Z'], star_i['omega_R'], star_i['omega_phi'], star_i['omega_z'], star_i['angle_R'], star_i['angle_phi'], star_i['angle_z'] = aAS.actionsFreqsAngles(
            #R,vR,vT,z,vz[,phi]
            star_i['R_Rzphi']*u.kpc,
            star_i['vR_Rzphi']*u.km/u.s,
            star_i['vT_Rzphi']*u.km/u.s,
            star_i['z_Rzphi']*u.kpc,
            star_i['vz_Rzphi']*u.km/u.s,
            star_i['phi_Rzphi']*u.rad,
            ro=r_galactic_centre,vo=v_circular
        )
    except:
        star_i['omega_R'] = [np.nan]
        star_i['omega_phi'] = [np.nan]
        star_i['omega_z'] = [np.nan]
        star_i['angle_R'] = [np.nan]
        star_i['angle_phi'] = [np.nan]
        star_i['angle_z'] = [np.nan]
        try:
            star_i['J_R'], star_i['L_Z'],star_i['J_Z'] = aAS(
                #R,vR,vT,z,vz[,phi]
                star_i['R_Rzphi']*u.kpc,
                star_i['vR_Rzphi']*u.km/u.s,
                star_i['vT_Rzphi']*u.km/u.s,
                star_i['z_Rzphi']*u.kpc,
                star_i['vz_Rzphi']*u.km/u.s,
                star_i['phi_Rzphi']*u.rad,
                ro=r_galactic_centre,vo=v_circular
            )
        except:
            star_i['J_R'] = [np.nan]
            star_i['L_Z'] = [np.nan]
            star_i['J_Z'] = [np.nan]

    try:
        star_i['ecc'], star_i['zmax'], star_i['R_peri'], star_i['R_ap'] = aAS.EccZmaxRperiRap(
            #R,vR,vT,z,vz[,phi]
            star_i['R_Rzphi']*u.kpc,
            star_i['vR_Rzphi']*u.km/u.s,
            star_i['vT_Rzphi']*u.km/u.s,
            star_i['z_Rzphi']*u.kpc,
            star_i['vz_Rzphi']*u.km/u.s,
            star_i['phi_Rzphi']*u.rad,
            ro=r_galactic_centre,vo=v_circular
        )         
        star_i['zmax']
        star_i['R_peri']
        star_i['R_peri']

    except:
        star_i['ecc'] = [np.nan]
        star_i['zmax'] = [np.nan]
        star_i['R_peri'] = [np.nan]
        star_i['R_ap'] = [np.nan]

    star_i['Energy'] = o.E(pot=pot,ro=r_galactic_centre,vo=v_circular,zo=z_galactic_plane)

    return(star_i)


# In[ ]:


def get_orbit_information_each_star(orbit_information, star_index):
    
    best_values = get_orbit_calculation_input(data = galah_gaia[star_index], MC_size = 1)
    star_i = estimate_orbit_parameters(star_index, orbit_calculation_input = best_values)

    for each_label in orbit_labels:
        try:
            orbit_information[each_label][star_index] = star_i[each_label][0]
        except:
            print('did not work for '+each_label)

    MC_values = get_orbit_calculation_input(data = galah_gaia[star_index], MC_size = MC_size)
    star_i = estimate_orbit_parameters(star_index, orbit_calculation_input = MC_values)    

    for each_label in orbit_labels:
        try:
            percentiles = np.percentile(star_i[each_label], q=[5,50,95])
            orbit_information[each_label+'_5'][star_index] = percentiles[0]
            orbit_information[each_label+'_50'][star_index] = percentiles[1]
            orbit_information[each_label+'_95'][star_index] = percentiles[2]
        except:
            print('did not work for '+each_label)
            
[get_orbit_information_each_star(orbit_information, star_index) for star_index in range(nr_stars)];


# # Compute orbit information

# In[ ]:


galah_dynamics = collections.OrderedDict()

if subset != -1:
    galah_dynamics['sobject_id'] = galah_gaia_input['sobject_id']
    galah_dynamics['use_dist_flag'] = galah_gaia_input['use_dist_flag']
    galah_dynamics['use_rv_flag'] = galah_gaia_input['use_rv_flag']
else:
    galah_dynamics['sobject_id'] = galah_gaia['sobject_id']
    galah_dynamics['use_dist_flag'] = galah_gaia['use_dist_flag']
    galah_dynamics['use_rv_flag'] = galah_gaia['use_rv_flag']

for each_orbit_label in orbit_labels:
    galah_dynamics[each_orbit_label] = np.zeros(nr_galah_stars, dtype=float)
    galah_dynamics[each_orbit_label].fill(np.nan)
    (galah_dynamics[each_orbit_label])[nr_galah_stars_dynamics] = orbit_information[each_orbit_label]

for each_orbit_label in orbit_labels:
    for each_sampler in ['5','50','95']:
        galah_dynamics[each_orbit_label+'_'+each_sampler] = np.zeros(nr_galah_stars, dtype=float)
        galah_dynamics[each_orbit_label+'_'+each_sampler].fill(np.nan)
        (galah_dynamics[each_orbit_label+'_'+each_sampler])[nr_galah_stars_dynamics] = orbit_information[each_orbit_label+'_'+each_sampler]

galah_dynamics_data = pandas.DataFrame(galah_dynamics,columns=galah_dynamics.keys())

data_for_fits = Table.from_pandas(galah_dynamics_data)
data_for_fits.write(out_dir+'dynamics_output/sobject_dynamics_'+str(subset)+'.fits',overwrite=True)


# In[ ]:




