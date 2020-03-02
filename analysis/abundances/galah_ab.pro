PRO GALAH_ab,field,object,setup,mode,cannon=cannon,ps=ps
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; This is the major program of the GALAH pipeline for step 2 (of 2):
; Elemental Abundance estimation
;
; It consists of 5 main parts: 
;
; 1) Initialise setup (variables, atmospheres etc.) 
; 2) Read in spectra + resolution (+ eventually correct telluric & skyline)
; 3) Read in stellar parameters
; 4) Optimisation of elemental abundances for
;      a) only one element with keyword mode='ZZ'
;      b) an element line with keyword mode='ZZ1234'
;      c) all elements with keyword mode='all'  
; 4.1) Segment selection, normalisation, synthesis with all elements
; 4.2) Synthesis only with element of choice (fixed segments and normalisation)
; 4.3) Iterative optimisation (fixed segments and normalisation)     
; 5) Clean up / diagnostic plots
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

; History
;
; 17/07/05 SB Starting from file used for reduction version 5.1
; 17/07/06 SB Adding comments to get better overview of (main) parts
; 17/07/06 SB Switch off telluric correction for IRAF reduction dr5.2
; 18/02/15 SB Set Li line.depth=0.99, because DEPTH files run /w sol. Li
; 18/02/28 SB Go to IRAF reduction 5.3 
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 1/5:
; Initialise setup
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@sme_struct

; Print and remember begin of calculations
start_time=systime(/SECONDS)
print, 'Starting GALAH_ab: ', systime()

; For virtual machine mode
cla=command_line_args(count=count)
if count ne 0 then begin
   field  = cla[0]
   object = cla[1]
   setup  = cla[2]
   mode   = cla[3]
endif

; Ensure strings without spaces 
field     = strcompress(field,/remove)
field_end = strmid(field,strlen(field)-4,4)
object    = strcompress(object,/remove)
setup     = strcompress(setup,/remove)
mode      = strcompress(mode,/remove)

if mode eq 'all' then begin
   readcol,'mode_'+setup,modeun,format='a',skipline=2,comment=';'
   modeun=unique(modeun)
endif else modeun=mode

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
obs_name       = field+'_'+object+'_'+setup+'_'+mode
obs_file       = obs_name+'.dat'
auto_alpha     = 0
depthmin       = 0.05d0 & print,'Using minimum depth ',depthmin
object_pivot   = object
NaN            = abs(sqrt(-1))
LBOL           = NaN
glob_free      = '-1'
cont_mask      = [setup+'_Cont.dat']                                                                                ;Dummy, not used 
segm_mask      = [obs_name+'_Segm.dat']                                                                             ;Small segments for each element
line_mask      = [obs_name+'.dat']                                                                                  ;First guess line mask
broad_lines    = [4861.3230d0,6562.7970d0, 5167.3216d,5172.6843d,5183.6042d]                                        ;Hbeta, Halpha, Mg Ib triplet
line_cores     = [8498.0200d0]                                                                                      ;example

elstr=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cs','Es']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Choose atmosphere and linelist
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
atmogrid_file  = 'marcs2014.sav'
   if atmogrid_file  eq 'marcs2014.sav'               then line_list      = 'galah_master_v5.2.fits'
   if atmogrid_file  eq 'marcs2012.sav'               then line_list      = 'galah_master_v5.2.fits'
   if atmogrid_file  eq 'stagger-t5havni_marcs.sav'   then line_list      = 'galah_master_v5.fits'
   if atmogrid_file  eq 'stagger-t5havni_marcs64.sav' then line_list      = 'galah_master_v5.fits'
   if atmogrid_file  eq 'stagger-tAA.sav'             then line_list      = 'galah_master_v5.fits'

print,'Using atmosphere grid: ',atmogrid_file
print,'Using linelist:        ',line_list


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 2/5:
; Read in spectra and possibly perform skyline and tellurics
; correction
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in 4 band spectra (ext. 0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
com='com' ; Assume the spectrum is not a multivisit/stacked spectrum 
if strmid(object,11,1) eq '2' then com='com2'
f1='SPECTRA/dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'1.fits'
f1_exists = file_test(f1)
f2='SPECTRA/dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'2.fits'
f3='SPECTRA/dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'3.fits'
f4='SPECTRA/dr5.3/'+strmid(object,0,6)+'/standard/'+com+'/'+object+'4.fits'

if not f1_exists then begin
   print,'Cannot find FITS file in SPECTRA/dr5.3: '+object+', will look in SPECTRA/'+field
   f1='SPECTRA/'+field+'/'+object+'1.fits'
   f1_exists = file_test(f1)
   f2='SPECTRA/'+field+'/'+object+'2.fits'
   f3='SPECTRA/'+field+'/'+object+'3.fits'
   f4='SPECTRA/'+field+'/'+object+'4.fits'

   if not f1_exists then begin
      print,'Cannot find FITS file in SPECTRA/'+field+': '+object
      goto,finishline2
   endif
endif

sob1  = mrdfits(f1,0,h1,/silent) & sob2  = mrdfits(f2,0,h2,/silent) & sob3  = mrdfits(f3,0,h3,/silent) & sob4  = mrdfits(f4,0,h4,/silent)
wave1 = fitswave(h1)             & wave2 = fitswave(h2)             & wave3 = fitswave(h3)             & wave4 = fitswave(h4)
serr1 = mrdfits(f1,1,h1,/silent) & serr2 = mrdfits(f2,1,h2,/silent) & serr3 = mrdfits(f3,1,h3,/silent) & serr4 = mrdfits(f4,1,h4,/silent)

wave  = [wave1, wave2, wave3, wave4]
sob   = [sob1 , sob2 , sob3 , sob4 ]
err   = [serr1, serr2, serr3, serr4]

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

iraf_dr53=mrdfits('DATA/sobject_iraf_53_2MASS_GaiaDR2_WISE_PanSTARRSDR1_BailerJones_K2seis_small.fits',1,/silent)
i = where(iraf_dr53.sobject_id eq long64(object),ic) & i = i[0]
if ic eq 0 then begin
   print,'Could not find entry in sobject_iraf53, return'
   goto,finishline2
endif else begin
   v_bary = iraf_dr53[i].v_bary
endelse
if i eq n_elements(iraf_dr53)-1 then print,'OBS: The correspnding entry in sobject_iraf*.fits is the last entry. Could not be found?'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Telluric correction (incr. errors)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;telluric                   =
;mrdfits('DATA/telluric_noao_21k.fits',1,telluric_h,/silent)
;wave_telluric              = telluric.wave/(1d0-v_bary/299792.458)
;telluric_interp            =
;interpol(telluric.flux,wave_telluric,wave)
;belowzero                  = where(telluric_interp lt 0.81)
;aboveone                   = where(telluric_interp gt 0.998)
;telluric_interp[belowzero] = 0.81
;telluric_interp[aboveone]  = 1.0
;err                        = err-1.+1./(telluric_interp*5.0-4.0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Skyline correction (incr. errors)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sky_mask   = mrdfits('DATA/Skyspectrum_161105.fits',1,/silent)
wave_sky   = sky_mask.wave/(1d0-v_bary/299792.458)
sky_interp = interpol(sky_mask.sky,wave_sky,wave)
;err        = err+sky_interp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Final uncertainties UOB
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
uob   = abs(sob * err)
wave1 = 0. & wave2 = 0. & wave3 = 0. & wave4 = 0.
sob1  = 0. & sob2  = 0. & sob3  = 0. & sob4  = 0.
serr1 = 0. & serr2 = 0. & serr3 = 0. & serr4 = 0.
i=where(finite(sob)) & wave  = wave[i] & sob   = sob[i] & uob   = uob[i]

printcol,'SPECTRA/'+obs_file,wave,sob,uob,smooth(sob,10,/NAN),format=specformat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; MAIN PART 3.1/5:
; Read-in stellar parameters
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 1st try Cannon
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

if keyword_set(cannon) then begin
   print,'Initial values according to CANNON SETUP'
   inp='GALAH_cannon_'+field+'.fits'
   if not file_test(inp) then inp='GALAH_'+field+'.fits'
   if not file_test(inp) then begin
      print,'Neither GALAH_cannon_*_.fits not GALAH_*.fits found'
      goto,finishline2
   endif
   inp      = mrdfits(inp,1,h,/silent)
   objects  = inp.object
   ind      = where(object eq objects,indc) & ind=ind[0]
   if indc eq -1 then begin
      print,'Not found in GALAH_cannon_'+field+'.fits'
      goto,finishline2
   endif
   teff     = inp[ind].teff
   grav     = inp[ind].logg
   feh      = inp[ind].feh
   vradglob = inp[ind].vel
   if mode eq 'Fe' then begin
      print,'Using Fe to get vmic+vsini from Cannon output'
      turb    = vmic_vmac(teff,grav,feh)
      vmic    = turb[0]
      vsini   = turb[1]
      vmac    = 0.0
      print,'Cannon_results: TEFF, LOGG, FE/H, vrad', teff, grav, feh, vradglob
      marcs_abund,abund,fehmod=feh
   endif else begin
      print,'Using SME Fe run output'
      restore,'OUTPUT/'+str_replace(obs_name,setup+'_'+mode,setup+'_Fe')+'_SME.out'
      teff = sme.teff
      grav = sme.grav
      feh  = sme.feh
      vmic = sme.vmic
      vmac = sme.vmac
      vsini= sme.vsini
      abund= sme.abund
   endelse
endif else begin
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; 2nd try SME output and collected FITS
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

   print,'Initial values according to SME SETUP'

   sme_sp = 'OUTPUT/'+field+'_'+object+'_'+setup+'_Sp_SME.out'
   ;if not file_test(sme_sp) and field ne 'K_lbol' and field ne 'HIP_lbol' and field ne 'test' then begin
   ;   print,'No old SME.out found, returning'
   ;   goto,finishline2
   ;endif
   inp='GALAH_'+field+'.fits'
   if not file_test(inp) then begin
      print,'No GALAH_'+field+'.fits found, returning'
      goto,finishline2
   endif
   inp      = mrdfits(inp,1,h,/silent)
   objects  = inp.object
   ind      = where(object eq objects,indc) & ind=ind[0]
   if indc eq -1 then begin
      print,'Not found in GALAH_'+field+'.fits'
      goto,finishline2
   endif
   vradglob = inp[ind].vel
   if not finite(vradglob) then goto,finishline2
   teff     = inp[ind].teff
   if not finite(teff) then goto,finishline2
   grav     = inp[ind].logg
   feh      = inp[ind].feh
   turb     = vmic_vmac(teff,grav,feh)
   vmic     = turb[0]
   vsini    = inp[ind].vsini
   vmac     = inp[ind].vmac
   ;if field ne 'K_lbol' then begin
   ;   restore,sme_sp
   ;   abund    = sme.abund
   ;endif else begin 
   marcs_abund,abund,fehmod=feh
   ;endelse

endelse

; Increase initial Li abund
abund[3-1] = abund[3-1] - feh + 1.25
; Increase N abund for giants (logg < 2)
; SB: Taken out, because this is causing inconsistent Chi2 behaviour
;if grav lt 2.0 then abund[7-1] = abund[7-1] + 0.4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 4) MAIN PART 4/5, ELEMENT ABUNDANCE OPTIMISATION
;    This can happen for
;    a) only one element with keyword mode='ZZ'
;    b) an element line with keyword mode='ZZ1234'
;    c) all elements with keyword mode='all'
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iter=[1,1,20]

for m=0,n_elements(modeun)-1 do begin 

   start_time_mode=systime(/SECONDS)
   print, 'Starting GALAH_ab MODE: ', systime()

   mode=modeun[m]

   print,'Running element: ',mode

   if strlen(mode) gt 2 then elem=strmid(mode,0,strlen(mode)-4) else elem=mode
   atom=where(elstr eq elem)    ;Chose element from list
   atom=atom[0]+1 

   readcol,'mode_'+setup,modeall,line0,line_st,line_en,segm_st,segm_en,ipres,format='a,d,d,d,f,f,f',comment=';',/silent
   i=where(modeall eq mode,ic)
   if ic eq 0 then begin 
      print,mode+' not found in mode_'+setup
      goto,finishline
   endif 

   line0 = line0[i] & line_st=line_st[i] & line_en=line_en[i] & segm_st=segm_st[i] & segm_en=segm_en[i] & ipres=ipres[i]
   segm_st=unique(segm_st)
   segm_en=unique(segm_en)

   ;This part is to get the resolution from the resolution map                                                                    
   for r=0,n_elements(ipres)-1 do begin
      line_mid=0.5*(line_st[r]+line_en[r])
      ipres[r]=interpol(resy,resx,line_mid)*resolution_factor
   endfor

   ;Check if any line is left in line-list
   line=mrdfits('LINELIST/'+line_list,1,/silent)
   line.lambda = double(string(line.lambda,format='(d16.4)'))
   line_list = obs_name+'.fits'    
   
   print,'Using trilinearly interpolated depth'
   interp_depth_grid,depth
   line.depth = [depth(100471:162356),depth(318432:435589),depth(652638:748075),depth(851743:884262)]
   ; Increase depth and adjust log_gf for Si7800
   si_7800 = where(line.lambda eq 7799.9957d and line.name[0] eq 'Si')
   line[si_7800].depth = 0.99 & line[si_7800].log_gf = -0.75
   ; Increase depth for Li6707
   li_6707 = where(line.lambda eq 6707.7635d and line.name[0] eq 'Li')
   line[li_6707].depth = 0.99

   ;readcol,'DATA/'+line_mask,line0,line_st,line_en,/silent,format='d,d,d',comment=';'
   jcall=lonarr(n_elements(line0))
   for i=0,n_elements(line0)-1 do begin
      depth_check = where(abs(line.lambda-line0[i]) lt 1e-5 and strtrim(line.name[0],2)+strtrim(line.name[1],2) eq elem,jc)
      j=where(abs(line.lambda-line0[i]) lt 1e-5 and strtrim(line.name[0],2)+strtrim(line.name[1],2) eq elem and line.depth gt depthmin,jc)
      jcall[i]=jc
   endfor
   if total(jcall) eq 0 then begin
      print,'Predicted line(s) too weak, exiting'
      goto,finishline
   endif  

   ;Print new line and segment mask with visible lines
   i=where(jcall ne 0)
   line0=line0[i] & line_st=line_st[i] & line_en=line_en[i]
   printcol,'DATA/'+line_mask,line0,line_st,line_en,format=maskformat
   printcol,'DATA/'+segm_mask,segm_st,segm_en,ipres,format=segmformat

   ;readcol,'DATA/'+segm_mask,segm_st,segm_en,format='d,d'
   for i=0,n_elements(line0)-1 do begin 
      j=where(segm_st lt line0[i] and segm_en gt line0[i]) & j=j[0]
      k=where(line.lambda gt segm_st[j] and line.lambda lt segm_en[j] and line.depth gt depthmin,jc)
      if i eq 0 then full=k else full=[full,k]
      k=where(line.lambda gt segm_st[j] and line.lambda lt segm_en[j] and strtrim(line.name[0],2)+strtrim(line.name[1],2) eq elem and line.depth gt depthmin,jc)
      if i eq 0 then clean=k else clean=[clean,k]
   endfor

   for i=0,n_elements(broad_lines)-1 do begin
      j=where(line.lambda eq broad_lines[i])
      clean=[clean,j]
      full=[full,j]
   endfor

   clean=unique(clean)
   full=unique(full)

   clean=clean[sort(clean)]
   full=full[sort(full)]

   mwrfits,line[full],'LINELIST/'+line_list,/create

   nlte=1 ; set to 0 for AVATAR runs
   nlte_elem_flags = bytarr(99)
   nlte_grids = strarr(99)
   nltee = ['Li','C','O','Na','Mg','Al','Si','K','Ca','Mn','Ba']
   nltez = [  3 , 6 , 8 , 11 , 12,  13 , 14 , 19, 20 , 25 , 56 ]
   inlte = where(nltee eq elem) & inlte=inlte[0]

   if nlte eq 1 and inlte ne -1 then begin
      
      print, 'NLTE on for '+nltee[inlte]
      nlte_elem_flags[nltez[inlte]-1]  = 1B
      nlte_grids[nltez[inlte]-1]       = ['Amarsi19_']+nltee[inlte]+['.grd']
      if nltee[inlte] eq 'Li' and teff lt 5750. then nlte_grids[nltez[inlte]-1] = ['Amarsi19_']+nltee[inlte]+['_t3800_5750.grd']
      if nltee[inlte] eq 'Li' and teff ge 5750. then nlte_grids[nltez[inlte]-1] = ['Amarsi19_']+nltee[inlte]+['_t5750_8000.grd']

   endif else begin

      print, 'LTE only'

   endelse

   for k=0,n_elements(iter)-1 do begin

      if k eq 0 then print,' STARTING LOOP 1 -> Normalization and full synthesis'
      if k eq 1 then print,' STARTING LOOP 2 -> Element synthesis'
      if k eq 2 then print,' STARTING LOOP 3 -> Element abundance optimisation (max. 20 iterations)'
      
      maxiter=iter[k]

      ; Print info for the run
      print,field+'_'+object,teff,grav,feh,vmac,vmic,vsini,vradglob,k,mode,format=infoformat
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;
      ; LOOP STEP 4.1: SEGMENT SELECTION AND NORMALISATION WITH FULL SYNTHESIS
      ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

      if k eq 0 then begin
         
         cscale_flag=1
         norm=1
         ab_free[atom-1]=0
         
         if mode eq 'H' and iter[k] le 1 then glob_free='-1'   else ab_free[atom-1]=0

         if version gt 4.00 then begin
            print,'using SME pipeline version '+string(version)
            make_struct,run=2,nlte=nlte,norm=norm,/newsme
         endif else begin
            print,'using SME pipeline version '+string(version)
            make_struct,run=2,nlte=nlte,norm=norm
         endelse

         if keyword_set(ps) then begin
            spawn,'cp -f OUTPUT/'+obs_name+'_SME.out '+obs_name+'_SME_full.out'
            inspect,field,object,setup,mode,/ps,/norm,yr=[-0.1,1.1],label=3
            spawn,'mv -f '+obs_name+'.ps '+obs_name+'_0.ps'    
         endif
            
      endif else begin
         
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ;
         ; LOOP STEP 4.2: SYNTHESIS WITH ONLY ELEMENT OF MODE
         ;
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         
         if k eq 1 then begin

            vradglob=0d0

            cscale_flag=-3
            norm=0
            ab_free[atom-1]=0
            
            ;Restore full synthesis from k=0
            restore,'OUTPUT/'+obs_name+'_SME.out'
            sme_full=sme
            sme_full_mob=sme_full.mob
            sme_full_wave=sme_full.wave
            sme_full_smod=sme_full.smod

            obs_file=obs_name+'.dat'
            printcol,'SPECTRA/'+obs_file,sme.wave,sme.sob,sme.uob,sme.smod,format=specformat
            
            ;Run with list only containing the element
            mwrfits,line[clean],'LINELIST/'+line_list,/create

            if version gt 4.00 then begin
               print,'using SME pipeline version '+string(version)
               make_struct,run=2,nlte=nlte,/newsme
            endif else begin
               print,'using SME pipeline version '+string(version)
               make_struct,run=2,nlte=nlte
            endelse

            if keyword_set(ps) then begin
               inspect,field,object,setup,mode,/ps,/norm,yr=[-0.1,1.1],label=3
               spawn,'mv -f '+obs_name+'.ps '+obs_name+'_1.ps' 
            endif
            
            ;Restore element synthesis
            restore,'OUTPUT/'+obs_name+'_SME.out'
            sme_clean=sme
            sme_clean_mob = sme_clean.mob
            sme_clean_wave = sme_clean.wave
            sme_clean_smod = sme_clean.smod

            if keyword_set(ps) then begin
               spawn,'cp -f OUTPUT/'+obs_name+'_SME.out '+obs_name+'_SME_clean.out'
            endif

            ;Double-check that the wavelength scales are identical
            if max(abs(sme_full.wave-sme_clean.wave)) gt 1d-4 then begin
               npix=min([n_elements(sme_full.wave),n_elements(sme_clean.wave)])
               sme_full.smod[0:npix-1]=interpol(sme_full.smod,sme_full.wave,sme_clean.wave[0:npix-1])
               sme_full.wave[0:npix-1]=sme_clean.wave[0:npix-1]
            endif

            ;Build up new line mask by comparing synthetic spectra
            dum1    = 0.d0  
            dum2    = 0.d0 
            dum3    = 0.d0
            chimax  = 5d-3
            eps     = 1d-2      ;Buffer
            
            for i=0,n_elements(line0)-1 do begin
               mob=where(sme_clean.wave ge line_st[i] and sme_clean.wave le line_en[i] and (sme_full.smod-sme_clean.smod)^2 lt chimax,mobc)
               print,'Clean MOB: '+string(mobc)+' with chimax '+string(chimax)
               print,(sme_full.smod[mob]-sme_clean.smod[mob])^2
               if mobc ge 4 then begin 
                  for j=0,mobc-1 do begin
                     if j eq 0 then begin
                        dum1=[dum1,sme_clean.wave[mob[j]]-eps]
                        dum2=[dum2,sme_clean.wave[mob[j]]+eps]
                        dum3=[dum3,line0[i]]
                     endif else begin
                        if mob[j] eq mob[j-1]+1 then begin
                           dum2(n_elements(dum2)-1)=sme_clean.wave[mob[j]]+eps 
                        endif else begin
                           dum1=[dum1,sme_clean.wave[mob[j]]-eps]
                           dum2=[dum2,sme_clean.wave[mob[j]]+eps]
                           dum3=[dum3,line0[i]]
                        endelse
                     endelse
                  endfor
               endif else begin
                  print,'Less than 5 MOB left: '
               endelse
            endfor

            ;No clean regions left
            if n_elements(dum3) le 1 then begin
               print,'No clean regions left'
               goto,finishline
            endif

            ;Print new line-mask for k=2 and revert to full list
            printcol,'DATA/'+obs_name+'.dat',dum3[1:*],dum1[1:*],dum2[1:*],format=maskformat
            mwrfits,line[full],'LINELIST/'+line_list,/create

         endif else begin
            
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;
            ; LOOP STEP 4.3: OPTIMISATION OF ABUNDANCE OF ELEMENT OF MODE
            ;
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            if mode eq 'Fe' and keyword_set(cannon) then begin
               glob_free=['VRAD','VSINI']
               ab_free[atom-1]=1
            endif
            if mode eq 'Ba5854' then glob_free=['VRAD']
            if mode eq 'Ca6718' or mode eq 'Ca5868' then glob_free=['VRAD']
            ;if mode eq 'Ba6497' then glob_free=['VRAD']
            if mode eq 'Si6722' then glob_free=['VRAD']
      
;            if mode eq 'Li' and object eq '140713000701340' then glob_free=['VRAD']
;            if mode eq 'Li6708' then glob_free=['VRAD']
            cscale_flag=-3
            norm=0
            ab_free[atom-1]=1
            
            if version gt 4.00 then begin
               print,'using SME pipeline version '+string(version)
               make_struct,run=1,nlte=nlte,/newsme
            endif else begin
               print,'using SME pipeline version '+string(version)
               make_struct,run=1,nlte=nlte
            endelse

            ;Add sme.full to access spectrum used for
            ;normalisation and the continuum masks
            cmsave,file='OUTPUT/'+obs_name+'_SME.out',sme_full_mob,sme_full_wave,sme_full_smod, sme_clean_mob, sme_clean_wave, sme_clean_smod,/append

            if keyword_set(ps) then begin
               inspect,field,object,setup,mode,/ps,/norm,yr=[-0.1,1.1],label=3
               spawn,'mv -f '+obs_name+'.ps '+obs_name+'_2.ps'
            endif
         endelse
      endelse
   endfor
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ; 4) MAIN PART 5/5, CLEAN UP AND DIAGNOSTICS
   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

   duration = floor(systime(/SECONDS) - start_time_mode)
   duration = strtrim( duration/3600,2)+'h' $
              +strtrim((duration mod 3600)/60,2)+'m' $
              +strtrim( duration mod 60,2)+'s'

   print, 'Finishing GALAH_ab after ',duration

   ;print,'Will inspect now'
   ;inspect,field,object,setup,mode,/ps,/norm,yr=[-0.1,1.1],label=3
   ;print,'Inspecting successful'
   ;spawn,'mv -f '+obs_name+'.ps OUTPUT'

   finishline:

   spawn,'rm -f DATA/'+obs_name+'.dat'
   spawn,'rm -f DATA/'+obs_name+'_Segm.dat'
   spawn,'rm -f LINELIST/'+obs_name+'.fits'
   spawn,'rm -f OUTPUT/'+obs_name+'_SME.inp'

   spawn,'rm -f SPECTRA/'+obs_file

   if m ne n_elements(modeun)-1 then print,'Continuing with next element/line'
endfor

finishline2:

duration = floor(systime(/SECONDS) - start_time)
duration = strtrim( duration/3600,2)+'h' $
           +strtrim((duration mod 3600)/60,2)+'m' $
           +strtrim( duration mod 60,2)+'s'
print, 'Finally Finishing GALAH_ab after ',duration

;If running in VM-mode, save the log-file
;spawn,'mv -f idl_'+obs_name+'.log OUTPUT/'
spawn,'mv -f idl_'+obs_name+'.log OUTPUT_'+field+'/'
spawn,'mv -f OUTPUT/'+obs_name+'_SME.out OUTPUT_'+field+'/'

END
