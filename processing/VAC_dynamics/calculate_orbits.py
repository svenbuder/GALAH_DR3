#!/usr/bin/env python
# coding: utf-8

# # Actions and Orbit caluclation with MC sampling for GALAH DR3
# 
# ## Author: Sven Buder
# 
# ### History:
# 181011 SB Created  
# 190222 SB Included sampling with 5D covariance matrix and fixed galpy coordinate transformation for J2015.5 in ICRS  
# 201001 SB Change to McMillan17 potential, including different RO and VO

# # What information you need
# 
# ra, dec, pmra, pmdec from Gaia DR2  
# 
# distance:  
# if you want to use parallax: parallax and parallax_uncertainty  
# if you want to use covariances: covariance entries from Gaia DR2  
# if you want to use Bailer-Jones distances: r_est, r_lo, r_hi  
# if you want to use BSTEP: dist_gbm, e_dist_gbm  
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

# Basic packages
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import os
import sys
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
from galpy.util import bovy_coords
from galpy.orbit import Orbit

# Reference values
#r_galactic_centre = 8.178*u.kpc # Gravity Collaboration, 2018, A&A, 615, 15
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


calculate_sun = False

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
    sun['vphi_Rzphi'] = o.vT()#*u.km/u.s        
    sun['vz_Rzphi'] = o.vz()#*u.km/u.s
    sun['vT_Rzphi'] = o.vT()#*u.km/u.s

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
    print('vR,vphi,vT,vz: '+"{:.2f}".format(sun['vR_Rzphi'])+' '+"{:.2f}".format(sun['vphi_Rzphi'])+' '+"{:.2f}".format(sun['vT_Rzphi'])+' '+"{:.2f}".format(sun['vz_Rzphi']))
    print('J_R,L_Z,J_Z: '+"{:.2f}".format(sun['J_R'][0])+' '+"{:.2f}".format(sun['L_Z'][0])+' '+"{:.2f}".format(sun['J_Z'][0]))
    print('ecc, zmax, R_peri, R_apo: '+"{:.2f}".format(sun['ecc'][0])+' '+"{:.2f}".format(sun['zmax'][0])+' '+"{:.2f}".format(sun['R_peri'][0])+' '+"{:.2f}".format(sun['R_ap'][0]))
    print('Energy: '+"{:.2f}".format(sun['Energy']))


# ### Input of 6D information in observable dimensions

# In[ ]:


try:
    galah_gaia_input = pyfits.getdata('/shared-storage/buder/svn-repos/trunk/GALAH/GALAH_DR3/catalogs/GALAH_DR3_main_200604_extended_caution_v2.fits',1)
except:
    galah_gaia_input = pyfits.getdata('/Users/svenbuder/GALAH_DR3/catalogs/GALAH_DR3_main_200604_extended_caution_v2.fits',1)

full_length = len(galah_gaia_input['sobject_id'])
print("Initial nr. of entries")
print(full_length)

subset_size = 10000

try:
    subset = int(sys.argv[1])
except:
    subset = 0
if subset*subset_size >= full_length:
    sys.exit('The subset is beyond the length of GALAH DR3')

galah_gaia_input = galah_gaia_input[subset*subset_size:np.min([(subset+1)*subset_size,full_length])]

nr_galah_stars = len(galah_gaia_input['sobject_id'])
print("Nr. stars per subset")
print(nr_galah_stars)

nr_galah_stars_dynamics = np.where(
    np.isfinite(galah_gaia_input['ra']) &
    np.isfinite(galah_gaia_input['dec']) &
    np.isfinite(galah_gaia_input['r_est']) &
    np.isfinite(galah_gaia_input['pmra']) &
    np.isfinite(galah_gaia_input['pmdec']) &
    #np.isfinite(galah_gaia_input['rv_guess']) &
    np.isfinite(galah_gaia_input['ra_error']) &
    np.isfinite(galah_gaia_input['dec_error']) &
    np.isfinite(galah_gaia_input['r_hi']) &
    np.isfinite(galah_gaia_input['r_lo']) &
    np.isfinite(galah_gaia_input['pmra_error']) &
    np.isfinite(galah_gaia_input['pmdec_error']) &
    #np.isfinite(galah_gaia_input['e_rv_guess']) &
    #(galah_gaia_input['rv_guess'] != 999.) &
    #(galah_gaia_input['rv_guess'] != 1024.) &
    (
        (
        np.isfinite(galah_gaia_input['rv_galah']) &
        np.isfinite(galah_gaia_input['e_rv_galah'])
        ) |
        (
        np.isfinite(galah_gaia_input['rv_gaia']) &
        np.isfinite(galah_gaia_input['e_rv_gaia'])
        )
    )
)[0]

# This should only me activated for tests with subsets of GALAH DR3
#nr_galah_stars_dynamics = nr_galah_stars_dynamics[:100]

galah_gaia = galah_gaia_input[nr_galah_stars_dynamics]
nr_stars = len(galah_gaia['sobject_id'])
print("Selected number of stars")
print(nr_stars)


# In[ ]:


six_dimensions = {}

# Right ascension [deg]
six_dimensions['ra'] = galah_gaia['ra']
# Declination [deg]
six_dimensions['dec'] = galah_gaia['dec']

# 1000./Parallax [mas]
six_dimensions['distance'] = 1000./galah_gaia['parallax']
# Bailer-Jones distance from Sun [pc]
six_dimensions['r_est'] = galah_gaia['r_est']
# BSTEP distance from Sun [pc]
six_dimensions['dist_gbm'] = galah_gaia['dist_gbm']*1000.
# Parallax [mas]
six_dimensions['parallax'] = galah_gaia['parallax']

# Total proper motion in direction of right ascension [mas/yr]
six_dimensions['pmra'] = galah_gaia['pmra']
# Total proper motion in direction of declination [mas/yr]
six_dimensions['pmdec'] = galah_gaia['pmdec']

# Radial velocity [km/s]
six_dimensions['vrad'] = galah_gaia['rv_galah']

# Use Gaia RVS if GALAH not good
use_gaia_instead = (
    (galah_gaia['e_rv_gaia'] < galah_gaia['e_rv_galah']) |
    (
        np.isnan(galah_gaia['e_rv_galah']) &
        np.isfinite(galah_gaia['e_rv_gaia'])
    )
)
six_dimensions['vrad'][use_gaia_instead] = galah_gaia['rv_gaia'][use_gaia_instead]


# In[ ]:


e_six_dimensions = {}

# Error of right ascension [mas] to [deg]
e_six_dimensions['ra'] = galah_gaia['ra_error']/(1000.*3600.)
# Error of declination [mas] to [deg]
e_six_dimensions['dec'] = galah_gaia['dec_error']/(1000.*3600.)
# Error of Bailer-Jones distance from Sun [pc]
e_six_dimensions['r_hi'] = galah_gaia['r_hi']
e_six_dimensions['r_lo']  = galah_gaia['r_lo']
    # We are currently sampling a 2-sided Gaussian because Bailer-Jones are only giving 16th/50th/86th percentiles.
    # Any idea how to improve that because of missing posteriors from Bailer-Jones?
# Error of BSTEP distance from Sun [pc]
e_six_dimensions['dist_gbm'] = galah_gaia['e_dist_gbm']*1000.
# Error of parallax [mas]
e_six_dimensions['parallax']  = galah_gaia['parallax_error']
# Error of total proper motion in direction of right ascension [mas/yr]
e_six_dimensions['pmra'] = galah_gaia['pmra_error']
# Error of total proper motion in direction of declination [mas/yr]
e_six_dimensions['pmdec'] = galah_gaia['pmdec_error']
# Error of radial velocity [km/s]
e_six_dimensions['vrad'] = galah_gaia['e_rv_galah']

# Use Gaia RVS if GALAH not good
use_gaia_instead = (
    (galah_gaia['e_rv_gaia'] < galah_gaia['e_rv_galah']) |
    (
        np.isnan(galah_gaia['e_rv_galah']) &
        np.isfinite(galah_gaia['e_rv_gaia'])
    )
)
e_six_dimensions['vrad'][use_gaia_instead] = galah_gaia['e_rv_gaia'][use_gaia_instead]


# ## Monte Carlo sampling of Orbits

# In[ ]:


MC_size = 10000

np.random.seed(123)

XYZ_labels       = ['X_XYZ','Y_XYZ','Z_XYZ']
UVW_labels       = ['U_UVW','V_UVW','W_UVW']

Rphiz_labels     = ['R_Rzphi','z_Rzphi','phi_Rzphi']
vRphiz_labels    = ['vR_Rzphi','vz_Rzphi','vphi_Rzphi','vT_Rzphi']

action_labels    = ['J_R','L_Z','J_Z']
ext_orbit_labels = ['ecc', 'zmax', 'R_peri', 'R_ap', 'Energy']

orbit_labels = np.concatenate((
    XYZ_labels,
    UVW_labels,
    Rphiz_labels,
    vRphiz_labels,
    action_labels,
    ext_orbit_labels    
    ))
print(orbit_labels)


# ### Samples

# In[ ]:


def sample_6d_uncertainty(
    six_dimensions,
    e_six_dimensions,
    MC_size=MC_size,
    use_BailerJones = False,
    use_BSTEP = False,
    parallax_offset=-0.029
    ):
    
    """
    This function samples the 6D space with the given uncertainties.
    4 Options are available:
    
    if MC_size==1: assume no uncertainties
    
    if use_BailerJones==True: Sample 6D parameters independently with distance from Bailer-Jones
    
    if no_correlation==False: Use Gaia DR2 covariance matrix to sample 5D
    and GALAH vrad for 6th D
    """

    np.random.seed(123)
    
    MC_sample_6D = {}
    
    # Option 1: We assume no errors and simply return the actual parameters
    if MC_size == 1:
        print('We assume no errors and simply return the actual parameters')
        for each_key in six_dimensions.keys():
            if each_key == 'distance':
                if use_BailerJones:
                    print('Using Bailer-Jones')
                    MC_sample_6D['distance'] = np.array([[six_dimensions['r_est'][x]] for x in range(nr_stars)])/1000.
                elif use_BSTEP:
                    print('Using BSTEP, otherwise Bailer-Jones')
                    MC_sample_6D['distance'] = np.array([[six_dimensions['r_est'][x]] for x in range(nr_stars)])/1000.
                    bstep = np.array([[six_dimensions['dist_gbm'][x]] for x in range(nr_stars)])/1000.
                    bstep_available = np.isfinite(bstep)
                    MC_sample_6D['distance'][bstep_available] = bstep[bstep_available]
                else:
                    print('Parallax')
                    MC_sample_6D['distance'] = np.array([[1000./(six_dimensions['parallax'][x]-parallax_offset)] for x in range(nr_stars)])/1000.
            else:
                MC_sample_6D[each_key] = np.array([[six_dimensions[each_key][x]] for x in range(nr_stars)])
    
    elif use_BailerJones:
        # Option 2: Sampling the distances from Bailer-Jones assuming 2 separate Gaussian distributions
        print('Sampling the distances from Bailer-Jones assuming 2 separate Gaussian distributions')
        distance_sigma_lo  = np.array([np.abs(np.random.normal(loc = 0., scale = six_dimensions['r_est'] -  e_six_dimensions['r_lo'])) for i in range(MC_size)])
        distance_sigma_hi  = np.array([np.abs(np.random.normal(loc = 0., scale = e_six_dimensions['r_hi'] - six_dimensions['r_est'])) for i in range(MC_size)])
        select_lo_hi = np.array([(np.random.uniform(0, 1, size=nr_stars) < 0.5).astype(float) for x in range(MC_size)])

        MC_sample_6D['ra']       = np.array([np.random.normal(loc=six_dimensions['ra'], scale=e_six_dimensions['ra']) for i in range(MC_size)]).T
        MC_sample_6D['dec']      = np.array([np.random.normal(loc=six_dimensions['dec'], scale=e_six_dimensions['dec']) for i in range(MC_size)]).T
        MC_sample_6D['distance'] = (six_dimensions['r_est'] + select_lo_hi*distance_sigma_hi - (1-select_lo_hi)*distance_sigma_lo).clip(min=0).T/1000.
        MC_sample_6D['pmra']     = np.array([np.random.normal(loc=six_dimensions['pmra'], scale=e_six_dimensions['pmra']) for i in range(MC_size)]).T
        MC_sample_6D['pmdec']    = np.array([np.random.normal(loc=six_dimensions['pmdec'], scale=e_six_dimensions['pmdec']) for i in range(MC_size)]).T
        MC_sample_6D['vrad']     = np.array([np.random.normal(loc=six_dimensions['vrad'], scale=e_six_dimensions['vrad']) for i in range(MC_size)]).T

    elif use_BSTEP:
        # Option 3: Using BSTEP GBM distances wherever possible (need useful stellar parameters, Bailer Jones otherwise)
        # Then check which values are not finite
        bstep_available = np.isfinite(six_dimensions['dist_gbm']) & np.isfinite(e_six_dimensions['dist_gbm'])
        nr_bstep = len(six_dimensions['dist_gbm'][bstep_available])

        MC_sample_6D['ra']       = np.array([np.random.normal(loc=six_dimensions['ra'], scale=e_six_dimensions['ra']) for i in range(MC_size)]).T
        MC_sample_6D['dec']      = np.array([np.random.normal(loc=six_dimensions['dec'], scale=e_six_dimensions['dec']) for i in range(MC_size)]).T
        
        # First fill everything with BSTEP
        print('Using BSTEP GBM distances (available for '+str(nr_bstep)+')')
        MC_sample_6D['distance'] = np.array([np.random.normal(loc=six_dimensions['dist_gbm'], scale=e_six_dimensions['dist_gbm']) for i in range(MC_size)]).T/1000.

        # Fill the ones without finite BSTEP with Bailer-Jones
        print('No parameters available for '+str(nr_stars-nr_bstep)+', using Bailer Jones for those')
        distance_sigma_lo  = np.array([np.abs(np.random.normal(loc = 0., scale = six_dimensions['r_est'][~bstep_available] -  e_six_dimensions['r_lo'][~bstep_available])) for i in range(MC_size)])
        distance_sigma_hi  = np.array([np.abs(np.random.normal(loc = 0., scale = e_six_dimensions['r_hi'][~bstep_available] - six_dimensions['r_est'][~bstep_available])) for i in range(MC_size)])
        select_lo_hi = np.array([(np.random.uniform(0, 1, size=np.shape(distance_sigma_lo)[1]) < 0.5).astype(float) for x in range(MC_size)])
        MC_sample_6D['distance'][~bstep_available,:] = (six_dimensions['r_est'][~bstep_available] + select_lo_hi*distance_sigma_hi - (1-select_lo_hi)*distance_sigma_lo).clip(min=0).T/1000.

        MC_sample_6D['pmra']     = np.array([np.random.normal(loc=six_dimensions['pmra'], scale=e_six_dimensions['pmra']) for i in range(MC_size)]).T
        MC_sample_6D['pmdec']    = np.array([np.random.normal(loc=six_dimensions['pmdec'], scale=e_six_dimensions['pmdec']) for i in range(MC_size)]).T
        MC_sample_6D['vrad']     = np.array([np.random.normal(loc=six_dimensions['vrad'], scale=e_six_dimensions['vrad']) for i in range(MC_size)]).T

    else:
        # Option4: We sample the errors including the covariance matrix
        print('We sample the errors including the covariance matrix and parallax offset')

        # Mean vector and covariance matrix
        mu = np.array(
            [six_dimensions['ra'],
             six_dimensions['dec'],
             six_dimensions['parallax']-parallax_offset,
             six_dimensions['pmra'],
             six_dimensions['pmdec'],
             six_dimensions['vrad']
            ])

        s00 = (e_six_dimensions['ra'])**2
        s11 = (e_six_dimensions['dec'])**2
        s22 = e_six_dimensions['parallax']**2
        s33 = e_six_dimensions['pmra']**2
        s44 = e_six_dimensions['pmdec']**2
        s55 = e_six_dimensions['vrad']**2

        s01 = (e_six_dimensions['ra']) * e_six_dimensions['dec'] * galah_gaia['ra_dec_corr']
        s02 = (e_six_dimensions['ra']) * e_six_dimensions['parallax'] * galah_gaia['ra_parallax_corr']
        s03 = (e_six_dimensions['ra']) * e_six_dimensions['pmra'] * galah_gaia['ra_pmra_corr']
        s04 = (e_six_dimensions['ra']) * e_six_dimensions['pmdec'] * galah_gaia['ra_pmdec_corr']
        s05 = 0.*np.ones(np.shape(galah_gaia['sobject_id'])[0])

        s12 = (e_six_dimensions['dec']) * e_six_dimensions['parallax'] * galah_gaia['dec_parallax_corr']
        s13 = (e_six_dimensions['dec']) * e_six_dimensions['pmra'] * galah_gaia['dec_pmra_corr']
        s14 = (e_six_dimensions['dec']) * e_six_dimensions['pmdec'] * galah_gaia['dec_pmdec_corr']
        s15 = 0.*np.ones(np.shape(galah_gaia['sobject_id'])[0])

        s23 = e_six_dimensions['parallax'] * e_six_dimensions['pmra'] * galah_gaia['parallax_pmra_corr']
        s24 = e_six_dimensions['parallax'] * e_six_dimensions['pmdec'] * galah_gaia['parallax_pmdec_corr']
        s25 = 0.*np.ones(np.shape(galah_gaia['sobject_id'])[0])

        s34 = e_six_dimensions['pmra'] * e_six_dimensions['pmdec'] * galah_gaia['pmra_pmdec_corr']
        s35 = 0.*np.ones(np.shape(galah_gaia['sobject_id'])[0])

        s45 = 0.*np.ones(np.shape(galah_gaia['sobject_id'])[0])

        sigma = np.array([
            [
            [s00[x], s01[x], s02[x], s03[x], s04[x], s05[x]],
            [s01[x], s11[x], s12[x], s13[x], s14[x], s15[x]],
            [s02[x], s12[x], s22[x], s23[x], s24[x], s25[x]],
            [s03[x], s13[x], s23[x], s33[x], s34[x], s35[x]],
            [s04[x], s14[x], s24[x], s34[x], s44[x], s45[x]],
            [s05[x], s15[x], s25[x], s35[x], s45[x], s55[x]]
            ] for x in range(np.shape(galah_gaia['sobject_id'])[0])
        ])

        sample = np.array([np.random.multivariate_normal(mu[:,x], sigma[x], size= MC_size) for x in range(np.shape(mu)[1])])

        print('Created MC_sample_6D with (Nr. entries, Nr. Samples, Dimensions):')
        print(np.shape(sample))

        MC_sample_6D['ra']       = sample[:,:,0] #in deg   #*np.pi/180. # in rad
        MC_sample_6D['dec']      = sample[:,:,1] #in deg   #*np.pi/180. # in rad
        MC_sample_6D['distance'] = 1./(sample[:,:,2]).clip(min=0.00001) # in kpc
        MC_sample_6D['pmra']     = sample[:,:,3] # in mas/yr
        MC_sample_6D['pmdec']    = sample[:,:,4] # in mas/yr
        MC_sample_6D['vrad']     = sample[:,:,5] # in km/s

    return MC_sample_6D


# # Compute orbit information

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


def estimate_orbit_parameters(MC_sample_6D, orbit_information, nr_stars):
    """
    Estimate orbit parameters from the given
    MC sample of 6D information for the Nr of stars
    and save it into orbit_information
    """

    for each_star in range(nr_stars):
        
        # We are creating a dictionary for each star
        star_i = dict()
        
        ra     = MC_sample_6D['ra'][each_star]           *u.deg
        dec    = MC_sample_6D['dec'][each_star]          *u.deg
        dist   = MC_sample_6D['distance'][each_star]     *u.kpc
        pm_ra  = MC_sample_6D['pmra'][each_star]         *u.mas/u.year
        pm_dec = MC_sample_6D['pmdec'][each_star]        *u.mas/u.year
        v_los  = MC_sample_6D['vrad'][each_star]         *u.km/u.s

        # Create the Orbit instance
        o = Orbit(
            vxvv=[ra,dec,dist,pm_ra, pm_dec,v_los],
            ro=r_galactic_centre,
            vo=v_circular,
            zo=z_galactic_plane,
            solarmotion=[-11.1, 15.17, 7.25]*u.km/u.s,
            radec=True
        )
        
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
        star_i['vphi_Rzphi'] = o.vT()#*u.km/u.s        
        star_i['vz_Rzphi'] = o.vz()#*u.km/u.s
        star_i['vT_Rzphi'] = o.vT()#*u.km/u.s
        
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

        for each_label in orbit_labels:
            if len(MC_sample_6D['ra'][0]) == 1:
                try:
                    orbit_information[each_label][each_star] = star_i[each_label][0]
                except:
                    print('did not work for '+each_label)
            else:
                try:
                    percentiles = np.percentile(star_i[each_label], q=[5,50,95])
                    orbit_information[each_label+'_5'][each_star] = percentiles[0]
                    orbit_information[each_label+'_50'][each_star] = percentiles[1]
                    orbit_information[each_label+'_95'][each_star] = percentiles[2]
                except:
                    print('did not work for '+each_label)
    return orbit_information


# In[ ]:


# We will first 'sample' only once with the best value
MC_sample_6D = sample_6d_uncertainty(six_dimensions,e_six_dimensions,MC_size=1,use_BSTEP=True)
orbit_information = estimate_orbit_parameters(MC_sample_6D, orbit_information, nr_stars)

# And now we sample with a certain Monte Carlo sampling size
MC_sample_6D = sample_6d_uncertainty(six_dimensions,e_six_dimensions,MC_size=MC_size,use_BSTEP=True)
orbit_information = estimate_orbit_parameters(MC_sample_6D, orbit_information, nr_stars)


# In[ ]:


# plot difference of distance estimation choices
if debug:
    MC_sample_6D = sample_6d_uncertainty(six_dimensions,e_six_dimensions,MC_size=MC_size,use_BSTEP=True)
    MC_sample_6D_1 = sample_6d_uncertainty(six_dimensions,e_six_dimensions,MC_size=MC_size,use_BailerJones=True)
    MC_sample_6D_2 = sample_6d_uncertainty(six_dimensions,e_six_dimensions,MC_size=MC_size)

    star_bla = 4

    print(galah_gaia['sobject_id'][star_bla],galah_gaia['teff'][star_bla],galah_gaia['logg'][star_bla],galah_gaia['fe_h'][star_bla],galah_gaia['flag_sp'][star_bla])

    if star_bla == 0:
        kwargs = dict(bins=np.linspace(0.16,0.1725,100),histtype='step')
    else:
        kwargs = dict(bins=50,histtype='step')

    print('Approach, best_dist, (best_dist-low_dist)/best_dist, low_dist, high_dist')
    print('BSTEP: ',six_dimensions['dist_gbm'][star_bla],e_six_dimensions['dist_gbm'][star_bla]/six_dimensions['dist_gbm'][star_bla],six_dimensions['dist_gbm'][star_bla]-e_six_dimensions['dist_gbm'][star_bla],six_dimensions['dist_gbm'][star_bla]+e_six_dimensions['dist_gbm'][star_bla])
    print('Bailer-Jones: ',six_dimensions['r_est'][star_bla],(six_dimensions['r_est'][star_bla]-e_six_dimensions['r_lo'][star_bla])/six_dimensions['r_est'][star_bla],e_six_dimensions['r_lo'][star_bla],e_six_dimensions['r_hi'][star_bla])
    print('Parallax: ',1000./six_dimensions['parallax'][star_bla],e_six_dimensions['parallax'][star_bla]/six_dimensions['parallax'][star_bla],1000./(six_dimensions['parallax'][star_bla]+e_six_dimensions['parallax'][star_bla]),1000./(six_dimensions['parallax'][star_bla]-e_six_dimensions['parallax'][star_bla]))
    print('BSTEP-Bailer-Jones',six_dimensions['dist_gbm'][star_bla]-six_dimensions['r_est'][star_bla],(six_dimensions['dist_gbm'][star_bla]-six_dimensions['r_est'][star_bla])/six_dimensions['dist_gbm'][star_bla])

    plt.hist(MC_sample_6D['distance'][star_bla],label='BSTEP',**kwargs);
    plt.hist(MC_sample_6D_1['distance'][star_bla],label='Bailer-Jones',**kwargs);
    plt.hist(MC_sample_6D_2['distance'][star_bla],label='Parallax',**kwargs);
    plt.legend()


# In[ ]:


def plot_sampling(data, star_index=1):
    f = plt.figure(figsize=(15,10))

    hist_kwarfs=dict(bins=25,cmin=1)
    hist_k = dict(bins=25)
    
    ax=plt.subplot(6,6,31)
    ax.hist(data['ra'][star_index],**hist_k);
    ax.set_xlabel('ra')
    ax.set_ylabel('ra')
    ax.set_xticks([])
    ax.set_yticks([])

    ax=plt.subplot(6,6,32)
    ax.hist2d(data['dec'][star_index],data['ra'][star_index],**hist_kwarfs);
    ax.set_xlabel('dec')
    ax.set_xticks([])
    ax=plt.subplot(6,6,33)
    ax.hist2d(data['distance'][star_index],data['ra'][star_index],**hist_kwarfs);
    ax.set_xlabel('distance')
    ax.set_yticks([])
    ax=plt.subplot(6,6,34)
    ax.hist2d(data['pmra'][star_index],data['ra'][star_index],**hist_kwarfs);
    ax.set_xlabel('pmra')
    ax.set_yticks([])
    ax=plt.subplot(6,6,35)
    ax.hist2d(data['pmdec'][star_index],data['ra'][star_index],**hist_kwarfs);
    ax.set_xlabel('pmdec')
    ax.set_yticks([])
    ax=plt.subplot(6,6,36)
    ax.hist2d(data['vrad'][star_index],data['ra'][star_index],**hist_kwarfs);
    ax.set_xlabel('rv')
    ax.set_yticks([])

    ax=plt.subplot(6,6,26)
    ax.set_ylabel('dec')
    ax.hist(data['dec'][star_index],**hist_k);
    ax.set_yticks([])
    ax=plt.subplot(6,6,27)
    ax.hist2d(data['distance'][star_index],data['dec'][star_index],**hist_kwarfs);
    ax.set_yticks([])
    ax=plt.subplot(6,6,28)
    ax.hist2d(data['pmra'][star_index],data['dec'][star_index],**hist_kwarfs);
    ax.set_yticks([])
    ax=plt.subplot(6,6,29)
    ax.hist2d(data['pmdec'][star_index],data['dec'][star_index],**hist_kwarfs);
    ax.set_yticks([])
    ax=plt.subplot(6,6,30)
    ax.hist2d(data['vrad'][star_index],data['dec'][star_index],**hist_kwarfs);
    ax.set_yticks([])

    ax=plt.subplot(6,6,21)
    ax.set_ylabel('parallax')
    ax.hist(data['distance'][star_index],**hist_k);
    ax=plt.subplot(6,6,22)
    ax.hist2d(data['pmra'][star_index],data['distance'][star_index],**hist_kwarfs);
    ax=plt.subplot(6,6,23)
    ax.hist2d(data['pmdec'][star_index],data['distance'][star_index],**hist_kwarfs);
    ax=plt.subplot(6,6,24)
    ax.hist2d(data['vrad'][star_index],data['distance'][star_index],**hist_kwarfs);

    ax=plt.subplot(6,6,16)
    ax.set_ylabel('pmra')
    ax.hist(data['pmra'][star_index],**hist_k);
    ax=plt.subplot(6,6,17)
    ax.hist2d(data['pmdec'][star_index],data['pmra'][star_index],**hist_kwarfs);
    ax=plt.subplot(6,6,18)
    ax.hist2d(data['vrad'][star_index],data['pmra'][star_index],**hist_kwarfs);

    ax=plt.subplot(6,6,11)
    ax.hist(data['pmdec'][star_index],**hist_k);
    ax.set_ylabel('pmdec')
    ax=plt.subplot(6,6,12)
    ax.hist2d(data['vrad'][star_index],data['pmdec'][star_index],**hist_kwarfs);

    ax=plt.subplot(6,6,6)
    ax.set_ylabel('vrad')
    ax.hist(data['vrad'][star_index],**hist_k);

    plt.tight_layout()

if debug==True:
    plot_sampling(MC_sample_6D, star_index = 0)


# In[ ]:


if debug==True:
    
    star_index = 2
    
    print("XYZ = ({x:8.2f},{y:8.2f},{z:8.2f}) [kpc]".format(
        x=orbit_information[XYZ_labels[0]][star_index],
        y=orbit_information[XYZ_labels[1]][star_index],
        z=orbit_information[XYZ_labels[2]][star_index]
        ))
          
    print("UVW = ({u:8.2f},{v:8.2f},{w:8.2f}) [kpc km/s]".format(
        u=orbit_information[UVW_labels[0]][star_index],
        v=orbit_information[UVW_labels[1]][star_index],
        w=orbit_information[UVW_labels[2]][star_index]
        ))

    print(r"R   = {r:6.2f} -{r_minus:6.2f} + {r_plus:6.2f} [kpc]".format(
            r=orbit_information[Rphiz_labels[0]][star_index],
            r_minus=orbit_information[Rphiz_labels[0]][star_index] - orbit_information[Rphiz_labels[0]+'_5'][star_index],
            r_plus=orbit_information[Rphiz_labels[0]+'_95'][star_index] - orbit_information[Rphiz_labels[0]][star_index]
            ))
    print(r"phi = {phi:6.2f} -{phi_minus:6.2f} + {phi_plus:6.2f} [kpc]".format(
            phi=orbit_information[Rphiz_labels[1]][star_index],
            phi_minus=orbit_information[Rphiz_labels[1]][star_index] - orbit_information[Rphiz_labels[1]+'_5'][star_index],
            phi_plus=orbit_information[Rphiz_labels[1]+'_95'][star_index] - orbit_information[Rphiz_labels[1]][star_index]
            ))
    print(r"z = {z:6.2f} -{z_minus:6.2f} + {z_plus:6.2f} [kpc]".format(
            z=orbit_information[Rphiz_labels[2]][star_index],
            z_minus=orbit_information[Rphiz_labels[2]][star_index] - orbit_information[Rphiz_labels[2]+'_5'][star_index],
            z_plus=orbit_information[Rphiz_labels[2]+'_95'][star_index] - orbit_information[Rphiz_labels[2]][star_index]
            ))

    print("J_R = {jr:6.2f} - {jr_m:6.2f} + {jr_p:6.2f} [kpc km/s]".format(
        jr=orbit_information['J_R'][star_index],
        jr_m=orbit_information['J_R'][star_index]-orbit_information['J_R_5'][star_index],
        jr_p=orbit_information['J_R_95'][star_index]-orbit_information['J_R'][star_index]
        ))
    print("L_Z = {lz:6.2f} - {lz_m:6.2f} + {lz_p:6.2f} [kpc km/s]".format(
        lz=orbit_information['L_Z'][star_index],
        lz_m=orbit_information['L_Z'][star_index]-orbit_information['L_Z_5'][star_index],
        lz_p=orbit_information['L_Z_95'][star_index]-orbit_information['L_Z'][star_index]
        ))
    print("J_Z = {jz:6.2f} - {jz_m:6.2f} + {jz_p:6.2f} [kpc km/s]".format(
        jz=orbit_information['J_Z'][star_index],
        jz_m=orbit_information['J_Z'][star_index]-orbit_information['J_Z_5'][star_index],
        jz_p=orbit_information['J_Z_95'][star_index]-orbit_information['J_Z'][star_index]
        ))
    print("e = {ecc:6.2f}, zmax = {zmax:6.2f}, Rperi = {rperi:6.2f}, Rapo = {rapo:6.2f}".format(
        ecc=orbit_information['ecc'][star_index],
        zmax=orbit_information['zmax'][star_index],
        rperi=orbit_information['R_peri'][star_index],
        rapo=orbit_information['R_ap'][star_index],
        ))


# In[ ]:


if debug==True: 
    
    useful = np.isfinite(galah_gaia['parallax'])
    #useful = (galah_gaia['parallax_error']/galah_gaia['parallax'] < 0.3) & (galah_gaia['parallax'] > 0)
    
    errorbar_kwargs = dict(fmt='o', rasterized=True, ms=1, c='r')
    
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3,figsize=(15,10))

    interim_orbit_information['Toomre'] = np.sqrt(interim_orbit_information['U_LSR']**2 + interim_orbit_information['W_LSR']**2)
    toomre_p = np.percentile(interim_orbit_information['Toomre'], q=[5,50,95], axis=1)
    
    orbit_information['Toomre'] = np.sqrt(orbit_information['U_LSR']**2 + orbit_information['W_LSR']**2)
    orbit_information['Toomre_5'] = toomre_p[0]
    orbit_information['Toomre_50'] = toomre_p[1]
    orbit_information['Toomre_95'] = toomre_p[2]
    
    def plot_distribution(xlabel, ylabel, yscale, ax):
        ax=ax
        if yscale == 'lin':
            ax.scatter(
                interim_orbit_information[xlabel][useful],
                interim_orbit_information[ylabel][useful],
                s=0.5,alpha=0.01,rasterized=True,
                label='MC'
            )
            ax.errorbar(
            orbit_information[xlabel+'_50'][useful],
            orbit_information[ylabel+'_50'][useful],
            xerr=[orbit_information[xlabel+'_50'][useful]-orbit_information[xlabel+'_5'][useful],orbit_information[xlabel+'_95'][useful]-orbit_information[xlabel+'_50'][useful]],
            yerr=[orbit_information[ylabel+'_50'][useful]-orbit_information[ylabel+'_5'][useful],orbit_information[ylabel+'_95'][useful]-orbit_information[ylabel+'_50'][useful]],
            label='5/50/95',
            **errorbar_kwargs
            )
            ax.scatter(
            orbit_information[xlabel][useful],
            orbit_information[ylabel][useful],
            c='k',rasterized=True,
            zorder=2,
            label='Best'
            )
        elif yscale == 'sqrt':
            ax.scatter(
                interim_orbit_information[xlabel][useful],
                np.sqrt(interim_orbit_information[ylabel][useful]),
                s=0.5,alpha=0.01,rasterized=True
            )
            ax.errorbar(
            orbit_information[xlabel+'_50'][useful],
            np.sqrt(orbit_information[ylabel+'_50'][useful]),
            xerr=[
                    orbit_information[xlabel+'_50'][useful]-orbit_information[xlabel+'_5'][useful],
                    orbit_information[xlabel+'_95'][useful]-orbit_information[xlabel+'_50'][useful]],
            yerr=[
                    np.sqrt(orbit_information[ylabel+'_50'][useful])-np.sqrt(orbit_information[ylabel+'_5'][useful]),
                    np.sqrt(orbit_information[ylabel+'_95'][useful])-np.sqrt(orbit_information[ylabel+'_50'][useful])
                ],
            **errorbar_kwargs
            )
            ax.scatter(
            orbit_information[xlabel][useful],
            np.sqrt(orbit_information[ylabel][useful]),
            c='k',rasterized=True,
            zorder=2
            )

    plot_distribution(xlabel=XYZ_labels[0], ylabel=XYZ_labels[1], yscale='lin', ax=ax1)
    plot_distribution(xlabel=Rphiz_labels[0], ylabel=Rphiz_labels[2], yscale='lin', ax=ax2)
    plot_distribution(xlabel=UVW_labels[0], ylabel='Toomre', yscale='lin', ax=ax3)
    plot_distribution(xlabel=ext_orbit_labels[0], ylabel=ext_orbit_labels[1], yscale='lin', ax=ax4)
    plot_distribution(xlabel=ext_orbit_labels[2], ylabel=ext_orbit_labels[3], yscale='lin', ax=ax5)
    plot_distribution(xlabel=action_labels[1], ylabel=action_labels[0], yscale='sqrt', ax=ax6)

    ax1.set_xlabel('X (XYZ) [kpc]')
    ax1.set_ylabel('Y (XYZ) [kpc]')

    ax2.set_xlabel('R (GC) [kpc]')
    ax2.set_ylabel('z (GC) [kpc]')
    legend = ax2.legend(loc='upper right',fontsize=15, markerscale=2)
    for ind, legend_handle in enumerate(legend.legendHandles):
        if ind==0:
            legend_handle.set_alpha(1)

    ax3.set_xlim(-700,200)
    ax3.set_ylim(0,450)
    ax3.set_xlabel('Toomre V (LSR) [km/s]')
    ax3.set_ylabel('Toomre UW (LSR) [km/s]')

    ax4.set_xlabel('Eccentricity')
    ax4.set_ylabel(r'$z_\text{max}$ [kpc]')

    ax5.set_xlabel(r'R (pericenter) [kpc]')
    ax5.set_ylabel(r'R (apocenter) [kpc]')

    ax6.set_xlabel(r'$L_Z$ [kpc km/s]')
    ax6.set_ylabel(r'$\sqrt{J_R \mathrm{[kpc km/s]}}$')
    ax6.set_ylim(-10,75)
    ax6.set_xlim(-1000,3000)

    plt.tight_layout()
    
    plt.savefig('figures/MC_output.pdf',dpi=300,bbox_inches='tight')


# In[ ]:


galah_dynamics = collections.OrderedDict()

galah_dynamics['sobject_id'] = galah_gaia_input['sobject_id']
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
data_for_fits.write('dynamics_output/sobject_dynamics_'+str(subset)+'.fits',overwrite=True)


# In[ ]:


try:
    import email
    import email.mime.application
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
    from email.MIMEImage import MIMEImage
    msg = MIMEMultipart()

    msg['From'] = 'gemini2'
    msg['To'] = 'buder@mpia.de'
    msg['Subject'] = 'Orbit calculation for subset '+str(subset)+' GALAH DR3 finished'                                                                                                         

    import smtplib
    mailer = smtplib.SMTP('localhost')
    mailer.sendmail('gemini2', 'buder@mpia.de', msg.as_string())
    mailer.close()
    print('Email sent')
except:
    print('Could not send email')


# # Posterior Sampling of the Distances/parallaxes

# In[ ]:


def likelihood_function(distance,parallax,uncertainty,parallax_zeropoint=-0.029/1000.):
    '''
    distance in pc
    parallax in mas
    uncertainty in mas
    
    returns unnormalized likelihood function
    '''
    parallax /= 1000
    uncertainty /= 1000
    result = np.exp(-np.divide((parallax-parallax_zeropoint-np.divide(1,distance))**2,2*uncertainty**2))
    return(result)

def prior_function(distance, r_len):
    prior = 1./(2.*r_len**3)*distance**2.*np.exp(-distance/r_len)
    prior[(distance<=0)] = 0.
    return prior

def posterior_function(prior, likelihood):
    return prior*likelihood

def cumulative_sum(x, y):
        """
        Compute areas of trapezoids and sum result
        """
        xdif = x[1:] - x[:-1]
        yavg = (y[0:-1] + y[1:] ) / 2.  
        tsum = np.cumsum(xdif*yavg) 

        return(tsum)

def calculate_posterior_percentiles(varpi, sigma_varpi, r_len, r_lo, r_est, r_hi):

    # sample within 5sigma parallax uncertainty
    distance_sample = 1000./np.linspace(
        varpi+5*sigma_varpi,
        (varpi-5*sigma_varpi).clip(min=0.1),
        1000)
    
    likelihood = likelihood_function(
        distance = distance_sample,
        parallax = varpi,
        uncertainty = sigma_varpi,
        parallax_zeropoint=-0.029/1000.)

    prior = prior_function(
        distance = distance_sample,
        r_len = r_len
    )
    
    posterior = posterior_function(prior, likelihood)
    
    cumulative_posterior = cumulative_sum(distance_linspace,prior*likelihood/np.max(prior*likelihood))
    cumulative_distance = (distance_linspace[1:]+distance_linspace[:-1])/2.
    cumulative_posterior /= np.max(cumulative_posterior)

    percentiles = np.interp([0.158655254,0.50,1-0.158655254],cumulative_posterior,cumulative_distance)

    plt.figure()
    plt.plot(
        distance_sample,
        posterior
    )
    for each in percentiles:
        plt.axvline(each,c='k')
    
    for each in [r_lo,r_est,r_hi]:
        plt.axvline(each,c='r')
        
    plt.show()
    
    return (percentiles)

# calculate_posterior_percentiles(
#     varpi = galah_gaia['parallax'][0], 
#     sigma_varpi = galah_gaia['parallax_error'][0], 
#     r_len = galah_gaia['r_len'][0], 
#     r_lo = galah_gaia['r_lo'][0], 
#     r_est = galah_gaia['r_est'][0], 
#     r_hi = galah_gaia['r_hi'][0]
# )

