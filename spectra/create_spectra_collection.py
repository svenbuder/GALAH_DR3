#!/usr/bin/env python
# coding: utf-8

# In[223]:


import numpy as np
import os
from os import path
import collections
from astropy.io import fits
from astropy.table import Table
from scipy.io.idl import readsav


# In[224]:


try:
    data = Table.read('DATA/GALAH_iDR3_main_test.fits')
except:
    data = Table.read('../catalogs/GALAH_iDR3_main_test.fits')


# In[227]:


try:
    mode = np.genfromtxt('mode_DR3',usecols=(0,),dtype=str,comments=';')
except:
    mode = np.genfromtxt('../input/mode_DR3',usecols=(0,),dtype=str,comments=';')
indexes = np.unique(mode, return_index=True)[1]
mode = [mode[index] for index in sorted(indexes)]


# In[228]:


def get_content(data_each, mode_each):
    """                                                                                                                                                                                                                   
    Get the spectrum information for each mode for each sobject_id into one                                                                                                                                               
    FITS file                                                                                                                                                                                                             
                                                                                                                                                                                                                          
    INPUT:                                                                                                                                                                                                                
    data_each : catalog entry from DR3 fits file                                                                                                                                                                          
    """
    content = collections.OrderedDict()
    content['wave'] = [np.NaN]
    content['sob'] = [np.NaN]
    content['uob'] = [np.NaN]
    content['smod'] = [np.NaN]

    file_exists = True

    try:
        sme = readsav((data_each['wg4_field']).replace(' ','')+'_'+str(data_each['sobject_id'])+'_DR3_'+mode_each+'_SME.out').sme[0]
        #sme = readsav('SMALL_OUTPUT/'+mode_each+'/OUTPUT_'+(data_each['wg4_field']).replace(' ','')+'/'+(data_each['wg4_field']).replace(' ','')+'_'+str(data_each['sobject_id'])+'_DR3_'+mode_each+'_SME.out').sme[0]
    except:
        try:
            sme = readsav('OUTPUT_'+(data_each['wg4_field']).replace(' ','')+'/'+(data_each['wg4_field']).replace(' ','')+'_'+str(data_each['sobject_id'])+'_DR3_'+mode_each+'_SME.out').sme[0]
        except:
            file_exists = False

    if file_exists:
        content['wave'] = sme.wave
        content['sob'] = sme.sob
        content['uob'] = sme.uob
        content['smod'] = sme.smod

    return(content)


# In[ ]:


for ind,data_each in enumerate(data):
    
    if ind%100==0:
        print(ind,len(data))

    new_hdul = fits.HDUList()

    cols = []
    
    for key in data_each.dtype.names:
        key_dtype = np.array(data_each[key]).dtype
        if str(key_dtype)[:2] == '<U':
            key_dtype = str(key_dtype).replace('<','')
        if str(key_dtype)[:1] == 'U':
            key_dtype = str(key_dtype).replace('U','A')
        else:
            key_dtype = key_dtype
            
        cols.append(fits.Column(name=key,format=key_dtype,array=[data_each[key]]))
        
    new_hdul.append(fits.BinTableHDU.from_columns(cols))

    for mode_each in mode:
        content = get_content(data_each, mode_each)

        col1 = fits.Column(name='wave',format='F',array=content['wave'])
        col2 = fits.Column(name='sob', format='F',array=content['sob'])
        col3 = fits.Column(name='uob', format='F',array=content['uob'])
        col4 = fits.Column(name='smod',format='F',array=content['smod'])

        new_hdul.append(fits.BinTableHDU.from_columns([col1,col2,col3,col4]))

    if not path.exists('SPECTRA_COLLECTION/'+str(data_each['sobject_id'])[:6]):
        os.system('mkdir SPECTRA_COLLECTION/'+str(data_each['sobject_id'])[:6])
        
    new_hdul.writeto('SPECTRA_COLLECTION/'+str(data_each['sobject_id'])[:6]+'/DR3_'+str(data_each['sobject_id'])+'.fits', overwrite=True)
    

