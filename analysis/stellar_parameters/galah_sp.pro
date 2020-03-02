PRO galah_sp,field,object,setup,cluster=cluster

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;
; This is the major program of the GALAH pipeline for step 1 (of 2):
; Stellar Parameter estimation 
;
; It consists of 5 main parts:
;
; 1) Initialise setup (variables, atmospheres etc.)
; 2) Read in spectra + resolution (+ eventually correct telluric & skyline)
; 3) Determine initial stellar parameters
; 3.1) Based on old Cannon/GUESS/default
; 3.2) If run with field_end LBOL, SEIS, FIXED, update initial parameters
; 3.3) Ensure reasonable parameters
; 4) Optimisation of stellar parameters (loop with alternating 4.1 and 4.2)
; 4.1) Segment selection and normalisation (fixed parameters)
; 4.2) Iterative optimisation (fixed segments and normalisation)
; 5) Clean up / diagnostic plots
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; History
;
; 17/07/05 SB Starting from file used for reduction version 5.1
; 17/07/06 SB Adding comments to get better overview of (main) parts
; 17/07/06 SB Switch off telluric correction for IRAF reduction dr5.2
; 18/02/28 SB Change to IRAF reduction dr5.3
; 18/02/28 SB Use marcs2014 atmo with new interp_atmo *pair and *grid
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 1/5:
; Initialise setup
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@sme_struct
  
; Print and remember begin of calculations
start_time=systime(/SECONDS)
print, 'Starting GALAH_sp: ', systime()
   
; For virtual machine mode   
cla=command_line_args(count=count)
if count ne 0 then begin
   field  = cla[0]
   object = cla[1]
   setup  = cla[2]
endif

; Ensure strings without spaces
field     = strcompress(field,/remove)
field_end = strmid(field,strlen(field)-4,4)
object    = strcompress(object,/remove)
setup     = strcompress(setup,/remove)

specformat     = '(d9.4,2x,3d25.8)'
infoformat     = '(a35,7f8.2,i4,a4)'
segmformat     = '(2d10.2,f10.0)'
maskformat     = '(d10.4,2d10.2)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize variables for later use
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
version        = 5.36
vrad_flag      = -1 ;-1: global vrad, 0: separate vrad for each wavelength segment  
cscale_flag    = 1 
vsini          = 2.0 
vradglob       = 0.0
ab_free        = intarr(99)
obs_name       = field+'_'+object+'_'+setup+'_Sp'
obs_file       = obs_name+'.dat'
auto_alpha     = 1
depthmin       = 0.1d0 & print,'Using minimum depth ',depthmin
broad_lines    = [4861.3230d0,6562.7970d0, 5167.3216d,5172.6843d,5183.6042d] ; Hbeta, Halpha, Mg Ib triplet
line_cores     = [8498.0200d0] ;example
object_pivot   = object
NaN            = abs(sqrt(-1))
LBOL           = NaN
b_mag = NaN & v_mag = NaN & g_mag = NaN & j_mag = NaN & h_mag = NaN & k_mag = NaN & w2_mag= NaN
plx   = NaN & dist  = NaN & e_dist= NaN & lbol  = NaN & e_lbol= NaN & ebv   = NaN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Choose atmosphere, grids and NLTE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
atmogrid_file  = 'marcs2014.sav'
   ; MARCS with Fe2016 grid and H2018 grid by Amarsi et al.
   if atmogrid_file  eq 'marcs2014.sav' then begin
      grid = ['marcs2012_H2018.grd','marcs2012_Fe2016.grd']
      line_list      = 'galah_master_v5.1.fits'
   endif
   ; STAGGER
;   if atmogrid_file  eq 'stagger-t5havni_marcs.sav' then begin
;      grid = 'NLTE/Fe_nlte_t5havni.grd'
;      line_list      = 'galah_master_v5.fits'
;   endif
;   if atmogrid_file  eq 'stagger-t5havni_marcs64.sav' then begin
;      grid = 'NLTE/Fe_nlte_t5havni64.grd'
;      line_list      = 'galah_master_v5.fits'
;   endif
;   if atmogrid_file  eq 'stagger-tAA.sav' then begin
;      grid = 'NLTE/Fe_nlte_t5havni.grd'
;      line_list      = 'galah_master_v5.fits'
;   endif

elem = [0,26]
nlte_elem_flags = bytarr(99) & nlte_elem_flags[elem-1] = 1B ; flag for NLTE computation
nlte_grids = strarr(99) & nlte_grids[elem-1] = grid         ; specify grids
nlte   = 1 ; Switches NLTE on (1) or off (0)

print,'Using atmosphere grid: ',atmogrid_file
print,'Using linelist:        ',line_list
print,'Using NLTE grid:       ',grid
if nlte eq 1 then print,'NLTE on' else print,'NLTE off'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 2/5:
; Read in spectra and possibly perform skyline and tellurics correction
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in 4 band spectra (ext. 0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
com='com' ; Assume the spectrum is not a multivisit/stacked spectrum
if strmid(object,11,1) eq '2' then com='com2'
if not keyword_set(cluster) then SPECTRA = 'SPECTRA/'
if keyword_set(cluster) then begin
   if cluster eq 'avatar' then SPECTRA = '/avatar/buder/trunk/GALAH/SPECTRA/'
   if cluster eq 'gemini2' then SPECTRA = '/shared-storage/buder/svn-repos/trunk/GALAH/SPECTRA'
endif
f1=SPECTRA+'dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'1.fits'
f1_exists = file_test(f1)
f2=SPECTRA+'dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'2.fits'
f3=SPECTRA+'dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'3.fits'
f4=SPECTRA+'dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'4.fits'

if not f1_exists then begin
   print,'Cannot find FITS file in SPECTRA/dr5.3: '+object+', will look in SPECTRA/'+field
   f1=SPECTRA+field+'/'+object+'1.fits'
   f1_exists = file_test(f1)
   f2=SPECTRA+field+'/'+object+'2.fits'
   f3=SPECTRA+field+'/'+object+'3.fits'
   f4=SPECTRA+field+'/'+object+'4.fits'

   if not f1_exists then begin 
      print,'Cannot find FITS file in SPECTRA/'+field+': '+object
      goto,finishline
   endif
endif

ir_exists = file_test(f4)

if ir_exists then begin
   sob1  = mrdfits(f1,0,h1,/silent) & sob2  = mrdfits(f2,0,h2,/silent) & sob3  = mrdfits(f3,0,h3,/silent) & sob4  = mrdfits(f4,0,h4,/silent)
   wave1 = fitswave(h1)             & wave2 = fitswave(h2)             & wave3 = fitswave(h3)             & wave4 = fitswave(h4)
   serr1 = mrdfits(f1,1,h1,/silent) & serr2 = mrdfits(f2,1,h2,/silent) & serr3 = mrdfits(f3,1,h3,/silent) & serr4 = mrdfits(f4,1,h4,/silent)

   wave  = [wave1, wave2, wave3, wave4]
   sob   = [sob1 , sob2 , sob3 , sob4 ]
   err   = [serr1, serr2, serr3, serr4]

endif else begin

   sob1  = mrdfits(f1,0,h1,/silent) & sob2  = mrdfits(f2,0,h2,/silent) & sob3   = mrdfits(f3,0,h3,/silent) ;& sob4  = mrdfits(f4,0,h4,/silent)
   wave1 = fitswave(h1)             & wave2 = fitswave(h2)             & wave3  = fitswave(h3)             ;& wave4 = fitswave(h4)
   serr1 = mrdfits(f1,1,h1,/silent) & serr2 = mrdfits(f2,1,h2,/silent) & serr3  = mrdfits(f3,1,h3,/silent) ;& serr4 = mrdfits(f4,1,h4,/silent)

   wave  = [wave1, wave2, wave3];, wave4]
   sob   = [sob1 , sob2 , sob3 ];, sob4 ]
   err   = [serr1, serr2, serr3];, serr4]

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; High-Res/Low-Res and Resolution Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
high_low = fxpar(h1,'SLITMASK')
if high_low eq 'IN      ' then begin
   resolution_factor = 1.789 ; increasing R from 28.000 to 50.000 for high-res targets
   print,'This spectrum is HIGH-RESOLUTION'
endif else begin
   resolution_factor = 1.000
   print,'This spectrum is LOW-RESOLUTION'
endelse

piv=fix(strmid(object,strlen(object)-3,strlen(object)))
res1=mrdfits('DATA/ccd1_piv.fits',0,res_h1,/silent)
res2=mrdfits('DATA/ccd2_piv.fits',0,res_h2,/silent)
res3=mrdfits('DATA/ccd3_piv.fits',0,res_h3,/silent)
res4=mrdfits('DATA/ccd4_piv.fits',0,res_h4,/silent)
resy1=res1[0:fxpar(res_h1,'NAXIS1')-1,piv-1] & resx1=fxpar(res_h1,'CRVAL1')+fxpar(res_h1,'CDELT1')*indgen(fxpar(res_h1,'NAXIS1'))
resy2=res2[0:fxpar(res_h2,'NAXIS1')-1,piv-1] & resx2=fxpar(res_h2,'CRVAL1')+fxpar(res_h2,'CDELT1')*indgen(fxpar(res_h2,'NAXIS1'))
resy3=res3[0:fxpar(res_h3,'NAXIS1')-1,piv-1] & resx3=fxpar(res_h3,'CRVAL1')+fxpar(res_h3,'CDELT1')*indgen(fxpar(res_h3,'NAXIS1'))
resy4=res4[0:fxpar(res_h4,'NAXIS1')-1,piv-1] & resx4=fxpar(res_h4,'CRVAL1')+fxpar(res_h4,'CDELT1')*indgen(fxpar(res_h4,'NAXIS1'))
resx=[resx1,resx2,resx3,resx4] & resy=[resy1,resy2,resy3,resy4]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in reduction pipeline output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;iraf_dr53=mrdfits('DATA/sobject_iraf_53_180227_vbary.fits',1,/silent)
;iraf_dr53=mrdfits('DATA/sobject_iraf_53_Gaia_2MASS_WISE.fits',1,/silent)
iraf_dr53=mrdfits('DATA/sobject_iraf_53_2MASS_GaiaDR2_WISE_PanSTARRSDR1_BailerJones_K2seis_small.fits',1,/silent)
in_dr53 = where(iraf_dr53.sobject_id eq long64(object),ic) & in_dr53 = in_dr53[0]
if ic eq 0 then begin
   print,'Could not find entry in sobject_iraf53, return'
   goto,finishline
endif else begin
   v_bary = iraf_dr53[in_dr53].v_bary   
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Telluric correction (incr. errors)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

telluric                   = mrdfits('DATA/telluric_noao_21k.fits',1,telluric_h,/silent)
wave_telluric              = telluric.wave/(1d0-v_bary/299792.458)
telluric_interp            = interpol(telluric.flux,wave_telluric,wave)
belowzero                  = where(telluric_interp lt 0.81)
aboveone                   = where(telluric_interp gt 0.998)
telluric_interp[belowzero] = 0.81
telluric_interp[aboveone]  = 1.0
err                        = err-1.+1./(telluric_interp*5.0-4.0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Skyline correction (incr. errors)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sky_mask   = mrdfits('DATA/Skyspectrum_161105.fits',1,/silent)
wave_sky   = sky_mask.wave/(1d0-v_bary/299792.458)
sky_interp = interpol(sky_mask.sky,wave_sky,wave)
err        = err+sky_interp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Final uncertainties UOB
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
uob   = sob * err
wave1 = 0. & wave2 = 0. & wave3 = 0. ;& wave4 = 0.
sob1  = 0. & sob2  = 0. & sob3  = 0. ;& sob4  = 0.
serr1 = 0. & serr2 = 0. & serr3 = 0. ;& serr4 = 0.
i=where(finite(sob)) & wave  = wave[i] & sob   = sob[i] & uob   = uob[i]

printcol,'SPECTRA/'+obs_file,wave,sob,uob,smooth(sob,10,/NAN),format=specformat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 3.1/5:
; Define initial parameters from Cannon/GUESS/default
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 1st try Cannon
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
cannon_version = 'iDR2'
cannon=mrdfits('DATA/sobject_iraf_iDR2_cannon_small.fits',1,/silent)
in_cannon=where(cannon.sobject_id eq object,ic) & in_cannon= in_cannon[0]

if ic eq 1 and cannon[in_cannon].flag_cannon le 1 then begin

   print,'Using '+cannon_version+' values as initial parameters'
   if cannon[in_cannon].flag_cannon ne 0 then begin
      print,'Continue, but Cannon flag is ',cannon[in_cannon].flag_cannon
   endif
   teff       = round(10.*cannon[in_cannon].teff_cannon)/10.
   grav       = round(1000.*cannon[in_cannon].logg_cannon)/1000.
   grav_input = grav
   feh        = round(1000.*cannon[in_cannon].feh_cannon)/1000.
   feh_input  = feh
   vmic       = cannon[in_cannon].vmic_cannon
   vsini      = cannon[in_cannon].vsini_cannon
   vmac       = 0.0
   vradglob   = round(100.*iraf_dr53[in_dr53].rv_guess)/100.
   print,field+'_'+object,teff,grav_input,feh_input,vradglob,'Input','Sp',format='(a35,4f8.2,a10,a4)'
   if teff lt 5250. and grav gt 4.25-(4.25-3.5)*(5250.-teff)/1250. then begin
      print,'Cool dwarf (upturn) identified, adjusting log(g) and [Fe/H]'
      grav = 4.5 + 0.2*(5250. - teff)/1250.
      feh  = 0
      turb=vmic_vmac(teff,grav,feh)
      vmic = turb[0]
      vsini = turb[1]
      vmac = 0.0
   endif
   ; Possible Parameters to distinguish between giants and dwarfs
   ; Estimated with Dartmouth isochrones for [Fe/H] = -0.5 dex, age = 15 Gyr 4500 K, 2.0 dex and 3650 K, 0.5 dex
   if teff lt 4500. and feh lt -0.75 and grav gt 2.0 - (2.0 - 0.5)*(4500.-teff)/850. then begin
      print,'Possible cool giant (TiO) identified, adjusting log(g) and [Fe/H]'
      ;grav = 4.5 + 0.2*(5250. - teff)/1250.
      feh  = 0
      turb=vmic_vmac(teff,grav,feh)
      vmic = turb[0]
      vsini = turb[1]
      vmac = 0.0
   endif
   if teff lt 4250 and grav lt 2.0 and feh lt -0.75 then begin
      print,'Giant at end of Cannon trainingset identified, adjusting [Fe/H]'
      feh = 0
   endif
   iter = [1,2,1,20] ;Four SME calls : normalise, iterate, normalise, iterate                                                                                                                  

endif else begin
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; 2nd try GUESS
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   print,'Bad initial values in '+cannon_version+' ('+string(ic)+' with flag '+string(cannon[in_cannon].flag_cannon)+'), trying GUESS values as initial parameters'
   if iraf_dr53[in_dr53].flag_guess eq 0 then begin
      teff=round(10.*iraf_dr53[in_dr53].teff_guess)/10.
      grav=round(1000.*iraf_dr53[in_dr53].logg_guess)/1000.
      grav_input=grav
      feh=round(1000.*iraf_dr53[in_dr53].feh_guess)/1000.
      feh_input=feh
      vradglob=round(100.*iraf_dr53[in_dr53].rv_guess)/100.
      print,field+'_'+object,teff,grav_input,feh_input,vradglob,'Input','Sp',format='(a35,4f8.2,a10,a4)'
      if grav gt 3.5 then begin
         grav=grav+0.5
         print,'NB: offset to initial GUESS log(g) by +0.5 for log(g)>3.5'
      endif
      feh=feh+0.15
      print,'NB: offset to initial GUESS [Fe/H] by 0.15'
      if teff lt 5250. and grav gt 4.25-(4.25-3.5)*(5250.-teff)/1250. then begin
         print,'Cool dwarf identified, adjusting log(g) and [Fe/H]'
         grav = 4.5 + 0.2*(5250. - teff)/1250.
         feh  = 0
      endif
      if teff lt 4250 and grav lt 2.0 and feh lt -0.75 then begin
         print,'Giant below 4250 K identified, adjusting [Fe/H] to 0.0'
         feh = 0
      endif
      turb=vmic_vmac(teff,grav,feh)
      vmic = turb[0]
      vsini = turb[1]
      vmac = 0.0
      iter = [1,2,1,20]         ;Four SME calls : normalise, iterate, normalise, iterate
      if iraf_dr53[in_dr53].red_flag ne 0 then begin
         print,'Reduction issue! ',iraf_dr53[in_dr53].red_flag ;bit mask for 1 to 8, summing um CCD problem, e.g. (0011) = 1+2=3 for BG
         if iraf_dr53[in_dr53].red_flag eq 14 then print,'Skyflat!'
      endif
   endif else begin
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; 3rd try default + Gaussian RV or return
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
      print,'Neither Cannon nor GUESS appropriate, trying Gaussian approach'
      rv_gaussfit,wave,sob,[4861.3230,6562.7970],range=[5,5],/auto,vradglob,e_vradglob
      if e_vradglob gt 5. then begin
         print,'Wether Cannon, GUESS nor Gaussian fit at Balmer lines are good, skipping to end'
         goto,finishline
      endif
      teff=5000.
      grav=3.0
      feh=-0.5
      turb = vmic_vmac(teff,grav,feh)
      vmic = turb[0]
      vsini = turb[1]
      vmac = 0.0
      print,field+'_'+object,teff,grav,feh,vradglob,'Input','Sp',format='(a35,4f8.2,a10,a4)'
      iter=[1,2,1,2,1,20]
   endelse
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 3.2/5:
; Account for special modes: LBOL, SEIS, FIXED
; If LBOL mode, prepare photometric/astrometric information
; If SEIS modee, prepare asteroseismic information
; If LBOL/SEIS mode, update LOGG based on the provided information
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
; If LBOL mode: GALAH+GaiaDR2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
if field_end eq 'lbol' then begin
   if iraf_dr53[in_dr53].r_est ge 0 then begin
      print,'Star in Gaia DR2 (and possibly 2MASS and WISE)'
      j_mag  = iraf_dr53[in_dr53].j_m    & e_j_mag  = iraf_dr53[in_dr53].j_msigcom
      h_mag  = iraf_dr53[in_dr53].h_m    & e_h_mag  = iraf_dr53[in_dr53].h_msigcom
      k_mag  = iraf_dr53[in_dr53].ks_m   & e_k_mag  = iraf_dr53[in_dr53].ks_msigcom
      w2_mag = iraf_dr53[in_dr53].w2mpro & e_w2_mag = iraf_dr53[in_dr53].w2mpro_error
      ebv = iraf_dr53[in_dr53].ebv
      a_k = 0.918*(h_mag - w2_mag - 0.08)
      if a_k lt 0.0 then a_k = 0.0
      print,'A_Ks from RJCE: '+string(a_k)+'+-'+string(sqrt(e_h_mag^2+e_w2_mag^2))
      print,'Quality flags (2MASS+WISE) '+string(iraf_dr53[in_dr53].ph_qual_tmass)+' & '+string(iraf_dr53[in_dr53].ph_qual_wise)
      print,'A_Ks from 0.36*E(B-V) '+string(0.36*ebv)
      if not finite(h_mag) or not finite(w2_mag) or strmid(string(iraf_dr53[in_dr53].ph_qual_tmass),1,1) ne 'A' or strmid(string(iraf_dr53[in_dr53].ph_qual_wise),1,1) ne 'A' then a_k=0.36*ebv
      if 0.36*ebv gt 3*a_k then begin
         ebv = 2.78*a_k
         print,'E(B-V) not trustworthy, setting to 2.78*A_Ks = '+string(ebv)+' instead of '+string(iraf_dr53[in_dr53].ebv)
      endif
      plx = iraf_dr53[in_dr53].parallax & e_plx = iraf_dr53[in_dr53].parallax_error
      dist = iraf_dr53[in_dr53].r_est
      print,'Distance   '+string(dist)+' with limits '+string(iraf_dr53[in_dr53].r_lo)+string(iraf_dr53[in_dr53].r_hi)
      print,'1/parallax '+string(1000./plx)+' and parallax: '+string(plx)+' with '+string(100*e_plx/plx)+' % uncertainty'
      e_dist = 0.5*(iraf_dr53[in_dr53].r_hi - iraf_dr53[in_dr53].r_lo)
      if dist lt 100. then begin
         print,'Star within 100pc, setting A_K and EBV to 0.0'
         a_k = 0.0
         ebv = 0.0
      endif
      if field eq 'Berkeley_73_lbol' then begin
         plx=0.128
         dist=6345.4
         e_dist=0.5*(17441.4-3882.9)
         e_plx = plx*(e_dist/dist)
         print,'New plx,e_plx,dist,e_dist by Cantat-Gaudin'
         print,plx,e_plx,dist,e_dist
      endif

   endif else begin
      print,'Star not in Gaia DR2, but lbol-keyword activated. Cancel.'
      goto,finishline
   endelse
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; if LBOL mode: Sun at 10 pc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
if field eq 'sun_lbol' or field eq 'RV_sun_lbol' then begin
   ; Here we assume the sun is at 10 pc, duhh!
   k_mag = 3.28 - 5.0 * ALOG10(1/10.) - 5.0
   h_mag = 3.32
   plx = 100.                   ; mas
   e_plx = 0.
   dist = 10. ;assuming sun is at 10 pc
   e_dist = 0.01
   ebv = 0.0
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; if SEIS mode: get numax
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if field_end eq 'seis' then begin
   i = where(iraf_dr53.sobject_id eq object and iraf_dr53.numax ge 0,ic) & i=i[0]
   if ic ne 0 then begin
      numax   = iraf_dr53[i].numax
      e_numax = 0. 
   endif else begin
      print,'No NUMAX value found, returning'
      goto,finishline
   endelse
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Update initial LOGG if LBOL/SEIS mode
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

if field_end eq 'lbol' then begin
   grav = update_logg(teff,grav,feh,/print) & grav=grav[0]
   turb = vmic_vmac(teff,grav,feh)
   vmic = turb[0] 
   vsini = turb[1] 
   vmac = 0.0
   print,'Running in LBOL mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, KMAG, PLX'
   print,'LBOL = '+string(lbol)+' Lbol_sun, KMAG = '+string(k_mag[0])+' mag, DISTANCE = '+string(dist[0])+', PLX = '+string(plx[0])+' mas, E_PLX/PLX = '+string(e_plx[0]/plx[0])
endif 
if field_end eq 'seis' then begin
   grav = update_logg(teff,grav,feh,/print)
   turb = vmic_vmac(teff,grav,feh)
   vmic = turb[0]
   vsini = turb[1]
   vmac = 0.0
   print,'Running in SEISMIC mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, NU_MAX'
   print,'NU_MAX = '+string(numax)+', E_NU_MAX/NU_MAX = '+string(e_numax/numax)
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Update TEFF/GRAV/FEH for FIXED mode
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if strpos(field,'fixed') ne -1 then begin
   readcol,'GALAH_'+field,fields,objects,setups,teffs,loggs,fehs,format=('a,a,a,f,f,f')
   i = where(objects eq object,ic) & i=i[0]
   teff = teffs[i]
   grav = loggs[i]
   feh  = fehs[i]
   print,'Running in FIXED mode, glob_free are VRAD, (VSINI)'
endif

if field_end ne 'lbol' and field_end ne 'seis' and strpos(field,'fixed') eq -1 then begin
   print,'Running in FREE mode, glob_free are TEFF, GRAV, FEH, VRAD, (VSINI)'
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 3.3/5:
; Ensure reasonable initial parameters
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Ensure reasonable parameters
teff=teff  < 7500. & teff=teff > 3000.
igrav=grav < 5.0   & grav=grav > 0.5
feh=feh    < 0.5   & feh=feh   > (-3.0)

; adjust teff and logg to marcs2012 lower grid limits
if atmogrid_file eq 'marcs2014.sav' then begin
   if teff lt 7500. and teff ge 6500. then grav=grav > 2.0
   if teff lt 6500. and teff ge 5750. then grav=grav > 1.5
   if teff lt 5750. and teff ge 5000. then grav=grav > 1.0
   if teff lt 5000. and teff ge 4000. then grav=grav > 0.5
   ; if teff lt 3000. then grav=grav > 0.5
endif

; adjust teff and logg to stagger grid limits
if atmogrid_file  eq 'stagger-t5havni_marcs64.sav' or atmogrid_file  eq 'stagger-t5havni_marcs.sav' or atmogrid_file eq 'stagger-tAA.sav' then begin
   if teff lt 4750. then grav=grav > 1.75
   if teff lt 5250. and teff ge 4750. then grav=grav > 2.75
   if teff lt 5750. and teff ge 5250. then grav=grav > 3.25
   if teff lt 6250. and teff ge 5750. then grav=grav > 3.75
   if teff ge 6250. then grav=grav > 4.15
endif

print,field+'_'+object,teff,grav,feh,vradglob,'Initial','Sp',format='(a35,4f8.2,a10,a4)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 4) MAIN PART 4/5, i.e. MAIN LOOP 'k': OPTIMISATION
; Even k / iter[k] eq 1: Normalisation
; Odd  k / iter[k] gt 1: Iterative optimisation
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

for k=0,n_elements(iter)-1 do begin 

   if k eq 0 then print,' STARTING LOOP 1 -> Pre-Normalization'
   if k eq 1 then print,' STARTING LOOP 2 -> First iterations (max. 2 iterations)'
   if k eq 2 then print,' STARTING LOOP 3 -> Normalization'
   if k eq 3 and n_elements(iter) eq 6 then print,' STARTING LOOP 4 -> Second iterations (max. 20 iterations)'
   if k eq 3 and n_elements(iter) eq 4 then print,' STARTING LOOP 4 -> Final iterations (max. 20 iterations)'
   if k eq 4 then print,' STARTING LOOP 5 -> Normalization'
   if k eq 5 then print,' STARTING LOOP 6 -> Final iterations (max. 20 iterations)'

   maxiter=iter[k]
   cont_mask   =  setup+'_Cont.dat'    ;Dummy continuum mask covering full segment
   segm_mask   =  setup+'_Segm.dat'    ;Segment definitions containing all lines

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   ; Initial abundances
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   if k eq 0 then begin 
      
      ;Initialise MARCS model abundances,
      ;N-enhanced in giants, Li assumed constant Spite Plateau, i.e. A(Li)=2.3
      marcs_abund,abund,fehmod=feh
      if grav lt 2.0 then abund[7-1] = abund[7-1] + 0.4
      ;abund[3-1] = abund[3-1] - feh + 1.25
   endif else begin

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Restore first loop and normalised spectrum
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      restore,'OUTPUT/'  + obs_name + '_SME.out'
      readcol,'SPECTRA/' + obs_file,wave,sob,uob,format='d,d,d',/silent

      ;Update parameters
      teff     = sme.teff 
      grav     = sme.grav
      if field_end eq 'seis' or field_end eq 'lbol' then grav = update_logg(teff,grav,feh,/print)
      feh      = sme.feh
      vradglob = sme.vrad
      abund    = sme.abund
      turb     = vmic_vmac(teff,grav,feh)
      vmic     = sme.vmic
      vsini    = sme.vsini
   endelse

   ;Print some information to screen
   print,field+'_'+object,teff,grav,feh,vmac,vmic,vsini,vradglob,k,'Sp',format=infoformat
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ; LOOP STEP 4.1: SEGMENT SELECTION AND NORMALISATION
   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   
   if iter[k] le 1 then begin
      ;No free parameters except normalisation using straight lines
      glob_free       = '-1'
      ab_free[*]      = 0
      cscale_flag     = 1

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Linelist, Depth-Grid, etc.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

      ;Read line-list and segments
      if k eq 0 then linefull = mrdfits('LINELIST/'+line_list,1)
      
      print,'Using trilinearly interpolated depth'
      interp_depth_grid,depth
      ;find_depth,depthfile
      ;depth=mrdfits(depthfile)
      if k eq 4 then print,'works'
      linefull.depth = [depth(100471:162356),depth(318432:435589),depth(652638:748075),depth(851743:884262)]
      
      if k eq 4 then print,'works1'
      line_list   =  obs_name+'.fits'
      line_mask   =  setup+'_Sp.dat' ;Line masks for log(g) sensitive lines, i.e. Fe,Ti,Sc
      
      readcol,'DATA/'+segm_mask,seg_st,seg_en,seg_ipres, comment=';',/silent,format='d,d,f'
      
      if k eq 4 then print,'works2'
      ;This part is to get the resolution from the resolution map
      for i=0,n_elements(seg_ipres)-1 do begin
         seg_mid=0.5*(seg_st[i]+seg_en[i])
         seg_ipres[i]=interpol(resy,resx,seg_mid)*resolution_factor
      endfor
      
      if k eq 4 then print,'works3'
      for i=0,n_elements(broad_lines)-1 do if i eq 0 then broad=where(linefull.lambda eq broad_lines[i]) else $
         broad=[broad,where(linefull.lambda eq broad_lines[i])]
      
      if k eq 4 then print,'works4'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      ;Loop over segments, normalise them, and update the line-list 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      for i=0,n_elements(seg_st)-1 do begin
         print,'Segment Nr. ',i+1,'/',n_elements(seg_st)
         
         ;Print a dummy fits-file, covering only one segment
         lines=where(linefull.lambda ge seg_st[i] and linefull.lambda le seg_en[i] and (linefull.depth gt depthmin or fs(linefull.name[0]) eq 'H'),lc)
         if lc eq 0 then lines=0
         if keyword_set(broad_lines) then plines=unique([lines,broad]) else plines=lines
         mwrfits,linefull[plines],'LINELIST/'+line_list,/create
         
         ;Print a dummy segment mask, covering only one segment
         printcol,'DATA/'+obs_name+'_Segm.dat',seg_st[i],seg_en[i],seg_ipres[i],format=segmformat
         segm_mask = obs_name+'_Segm.dat'
         
         ;Run SME for the segment and restore the output
         if k eq 0 then norm=1 else norm=0 ; pre-normalise for k eq 0
         
         if version ge 5.00 then begin
            print,version
            make_struct,run=2,norm=norm,nlte=nlte,/newsme
         endif else begin
            make_struct,run=2,norm=norm,nlte=nlte
         endelse

         restore,'OUTPUT/'+obs_name+'_SME.out'

         ;This should only be used when checking that the synthesis works
         ;plot,sme.wave,sme.sob,Color=cgColor('black'), Background=cgColor('white'),yrange=[-0.1,1.1]
         ;oplot,sme.wave,sme.smod,Color=cgColor('red')
         ;return
         
         ;Create dummy vectors to collect line-lists and normalised spectra
         if i eq 0 then begin
            linesave = plines
            dum1     = sme.wave
            dum2     = sme.sob
            dum3     = sme.uob
            dum4     = sme.smod
            dum5     = sme.mob
         endif else begin
            linesave = [linesave, plines ]
            dum1     = [dum1     ,sme.wave]
            dum2     = [dum2     ,sme.sob ]
            dum3     = [dum3     ,sme.uob ]
            dum4     = [dum4     ,sme.smod]
            dum5     = [dum5     ,sme.mob]
         endelse
         linesave = unique(linesave)
         
         ;Create special masks for Balmer lines
         balmer0 = [4861.3230d0,6562.7970d0]
         for b = 0,1 do begin
            if seg_st[i] lt balmer0[b] and seg_en[i] gt balmer0[b] then begin
               smef = sme
               line_list = 'galah_H.fits'
               print,'start make_struct Balmer'
               if version ge 5.00 then begin
                  print,version
                  make_struct,run=2,norm=norm,nlte=nlte,/newsme
               endif else begin
                  make_struct,run=2,norm=norm,nlte=nlte
               endelse

               print,'finished make_struct Balmer'
               restore,'OUTPUT/'+obs_name+'_SME.out'
               b0   = balmer0[b]
               balmer_mask,smef,sme,b0,balmer_st,balmer_en,teff,grav,feh
               line_list = obs_name+'.fits'
               readcol,'DATA/'+line_mask,line0,line_st,line_en,format='d,d,d',/silent
               line_mask = obs_name+'.dat' ;Line mask with Fe and blend-free H
               printcol,'DATA/'+line_mask,[line0,b0],[line_st,balmer_st],[line_en,balmer_en],format=maskformat
            endif
         endfor
      endfor
      
      ;Print full updated line-list
      mwrfits,linefull[linesave],'LINELIST/'+obs_name+'.fits',/create
      ;Print full normalised spectrum
      printcol,'SPECTRA/'+obs_file,dum1,dum2,dum3,dum4,format=specformat
      sme_full_smod=dum4
      sme_full_wave=dum1
      sme_full_mob=dum5
      dum1=0.
      dum2=0.
      dum3=0.
      dum4=0.

   ; END OF LOOP STEP 4.1
   endif else begin

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ; LOOP STEP 4.2: ITERATIVE OPTIMISATION
   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
   
      ; Free parameters
      if field_end eq 'lbol' then begin
         glob_free = ['TEFF', 'FEH', 'VRAD']
         print,'Running in LBOL mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, KMAG, PLX'
         print,'LBOL = '+string(lbol)+' Lbol_sun, KMAG = '+string(k_mag[0])+' mag, PLX = '+string(plx[0])+' mas, E_PLX/PLX = '+string(e_plx[0]/plx[0])
      endif else begin
         if field_end eq 'seis' then begin
            glob_free = ['TEFF', 'FEH', 'VRAD']
            print,'Running in SEISMIC mode, glob_free are TEF, FEH, VRAD, (VSINI), GRAV adjusted from TEFF, NU_MAX'
            print,'NU_MAX = '+string(numax)+', E_NU_MAX/NU_MAX = '+string(numax/e_numax)
         endif else begin
            if strpos(field,'fixed') ne -1 then begin
               glob_free = ['VRAD']
               print,'Running in FIXED mode, glob_free are VRAD, (VSINI)'
            endif else begin
               glob_free = ['TEFF', 'GRAV', 'FEH', 'VRAD']
               ;We could fit Sc, Ti, and Fe additionally 
               ;ab_free[20,21,25] = 1
               print,'Running in FREE mode, glob_free are TEFF, GRAV, FEH, VRAD, (VSINI)'
            endelse
         endelse
      endelse

      if vsini gt 1. then glob_free=[glob_free,'VSINI'] else vsini=1.
      
      line_list  = field+'_'+object+'_'+setup+'_Sp.fits'
      cscale_flag = -3   ;Normalisation fixed when solving for stellar parameters

     ;Initialise call to SME
     if version ge 5.00 then begin
        print,version
        make_struct,run=1,nlte=nlte,/newsme
     endif else begin
        make_struct,run=1,nlte=nlte
     endelse

     restore,'OUTPUT/'  + obs_name + '_SME.out'

   ; END OF LOOP STEP 4.2
   endelse

; END OF MAIN LOOP / MAIN PART 4
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 5/5:
; Clean up and diagnostic plot ('INSPECT')
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cmsave,file='OUTPUT/'+obs_name+'_SME.out',sme_full_mob,sme_full_wave,sme_full_smod,/append

duration = floor(systime(/SECONDS) - start_time)
  duration = strtrim( duration/3600,2)+'h' $
            +strtrim((duration mod 3600)/60,2)+'m' $
            +strtrim( duration mod 60,2)+'s'

print, 'Finishing GALAH_sp after ',duration

;Inspect SME result
;inspect,field,object,setup,'Sp',/ps,yr=[-0.1,1.1],label=5
;spawn,'mv -f '+obs_name+'.ps OUTPUT'

finishline:

;Clean up                                                                                                                                                                                                   

;spawn,'mv -f idl_'+obs_name+'.log OUTPUT/'
spawn,'rm -f DATA/'+obs_name+'.dat'
spawn,'rm -f SPECTRA/'+obs_name+'.dat'
spawn,'rm -f DATA/'+obs_name+'_Segm.dat'
spawn,'rm -f LINELIST/'+obs_name+'.fits'
spawn,'rm -f OUTPUT/'+obs_name+'_SME.inp'
test_dir = file_test('OUTPUT_'+field)
if not test_dir then spawn,'mkdir OUTPUT_'+field
spawn,'mv -f OUTPUT/'+obs_name+'_SME.out OUTPUT_'+field+'/'+obs_name+'_SME.out'
spawn,'mv -f idl_'+obs_name+'.log OUTPUT_'+field+'/'
;spawn,'rm -f core*'   

END
