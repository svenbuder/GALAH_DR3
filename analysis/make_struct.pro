PRO make_obs,status

;Reads observed spectrum and segment mask to create
;wave, sob, uob, wran, and nseg

@sme_struct
c=2.99792458d8 ; speed of light in m/s

status=0

;Read observations and synthesis
readcol,'SPECTRA/'+obs_file,dum1,dum2,dum3,dum5,/silent,format='d,d,d,d'

;Radial velocity correct 
dum1=dum1/(vradglob[0]*1d3/c + 1d0)

;Read segments and sort in increasing wavelength
readcol,'DATA/'+segm_mask,seg_st,seg_en,seg_ipres,comment=';',/silent,format='d,d,f'

;This sorts out segments 35
i=sort(seg_st) & seg_st=seg_st[i] & seg_en=seg_en[i] & seg_ipres=seg_ipres[i]

;Adjust resolution to resolution maps
piv=fix(strmid(object_pivot,strlen(object_pivot)-3,strlen(object_pivot)))
res1=mrdfits('DATA/ccd1_piv.fits',0,res_h1,/silent)
res2=mrdfits('DATA/ccd2_piv.fits',0,res_h2,/silent)
res3=mrdfits('DATA/ccd3_piv.fits',0,res_h3,/silent)
res4=mrdfits('DATA/ccd4_piv.fits',0,res_h4,/silent)
resy1=res1[0:fxpar(res_h1,'NAXIS1')-1,piv-1]
resx1=fxpar(res_h1,'CRVAL1')+fxpar(res_h1,'CDELT1')*indgen(fxpar(res_h1,'NAXIS1'))
resy2=res2[0:fxpar(res_h2,'NAXIS1')-1,piv-1]
resx2=fxpar(res_h2,'CRVAL1')+fxpar(res_h2,'CDELT1')*indgen(fxpar(res_h2,'NAXIS1'))
resy3=res3[0:fxpar(res_h3,'NAXIS1')-1,piv-1]
resx3=fxpar(res_h3,'CRVAL1')+fxpar(res_h3,'CDELT1')*indgen(fxpar(res_h3,'NAXIS1'))
resy4=res4[0:fxpar(res_h4,'NAXIS1')-1,piv-1]
resx4=fxpar(res_h4,'CRVAL1')+fxpar(res_h4,'CDELT1')*indgen(fxpar(res_h4,'NAXIS1'))
resx=[resx1,resx2,resx3,resx4]
resy=[resy1,resy2,resy3,resy4]    

;This part is to get the resolution from the resolution map
for i=0,n_elements(seg_ipres)-1 do begin
    seg_mid=0.5*(seg_st[i]+seg_en[i])
    seg_ipres[i]=interpol(resy,resx,seg_mid)*resolution_factor
endfor
print,'Resolution: ',seg_ipres

if min(seg_en-seg_st) le 0 then begin
    print,'Segment has negative range'
    return
endif

if n_elements(seg_en) gt 1 then begin
    if max(seg_en[0:*]-seg_st[1:*]) gt 0 then begin
        print,'Overlapping segments'
        return
    endif
endif

;If not normalising, make sure there are lines in the segment
k=0
if maxiter gt 1 then begin 
    readcol,'DATA/'+line_mask,line0,comment=';',/silent,format='d'
    for i=0,n_elements(seg_st)-1 do begin
        j=where(line0 gt seg_st[i] and line0 le seg_en[i],jc)
        if jc ne 0 then k=[k,i]
    endfor
    if n_elements(k) eq 1 then begin
        print,'No segments with lines'
        return
    endif
    seg_st=seg_st[k[1:*]]
    seg_en=seg_en[k[1:*]]
endif

;Select wavelength regions inside segments
k=0
for i=0,n_elements(seg_st)-1 do begin

   j=where(dum1 ge seg_st[i] and dum1 le seg_en[i], jc)
   if jc gt 10 then begin
       if k eq 0 then begin
           k = 1
           wave  =  dum1[j]
           sob   =  dum2[j]
           uob   =  dum3[j]
           ipres =  seg_ipres[i]
           smod  =  dum5[j]
           wind  =  [long(n_elements(wave)-1)] ; index of last pixel
       endif else begin
           wave  =  [wave  , dum1[j]]
           sob   =  [sob   , dum2[j]]
           uob   =  [uob   , dum3[j]]
           ipres =  [ipres , seg_ipres[i]]
           smod  =  [smod  , dum5[j]]
           wind  =  [wind  , long(n_elements(wave)-1)]
       endelse
   endif
endfor

if k eq 0 then begin
    print,'No observations in segment mask'
    return
endif

nseg = n_elements(wind)

;Define wran with exact segment ranges
wran = dblarr(2,nseg)

wran[0,0] = [wave[0]] 
wran[1,0] = [wave[wind[0]]]

for i=0,nseg-2 do begin
   wran[0,i+1]=[wave[wind[i]+1]] 
   wran[1,i+1]=[wave[wind[i+1]]]
endfor

status=1

END

PRO pre_norm

@sme_struct

;Performs first-guess normalisation by robustly converging straight line fit to high pixels

print,'Pre-normalising all segments'

for i = 0,nseg-1 do begin 

    j  = where(wave ge wran[0,i] and wave le wran[1,i] and sob gt 0.,jc)

    if jc gt 20 then begin
       co = autonorm(wave[j],sob[j],o=1,c=cont)
       sob[j] = sob[j]/cont
       uob[j] = uob[j]/cont
    endif else begin
       uob[j] = uob[j]/mean(sob[j])
       sob[j] = sob[j]/mean(sob[j])
    endelse

endfor

END

PRO make_mob

;Reads the line and continuum masks and defines mob

@sme_struct

mob = intarr(n_elements(wave))

; Read in mob data
readcol,'DATA/'+line_mask,line0,line_st,line_en,comment=';',/silent,format='d,d,d'
readcol,'DATA/'+cont_mask,cont_st,cont_en,comment=';',/silent,format='d,d'

; line mob - disregard regions with negative fluxes and above 1.1 or 1+3*sigma
;            At least one point must be <1.0

for i = 0,n_elements(line_st)-1 do begin
    m = where(wave ge line_st[i] and wave le line_en[i],mc)
    running_snr = []
    for ii =0,n_elements(sob)-1 do begin
       running_snr = [running_snr,max([1.0+10*1./mean(sob[max([0,ii-4]):min([ii+4,n_elements(sob)-1])]/uob[max([0,ii-4]):min([ii+4,n_elements(sob)-1])]),1.5])]
    endfor
    if mc ne 0 then begin 
;       if min(sob[m]) gt 0.0 and mean(sob[m]) le 1.0 and max(sob[m]) lt
;       max(running_snr[m]) then $ KL 06.08.2018
       if min(sob[m]) gt 0.0 and max(sob[m]) lt max(running_snr[m]) then $
          mob [m] = 1            
    endif
endfor

; cont mob  - continuum points selected between 0 and 1.2 or 1+3*sigma, where there are no line masks
;             avoid buffer zone at edges

if cscale_flag ge 0 then begin
    for i = 0,n_elements(cont_st)-1 do begin
        running_snr = []
        for ii =0,n_elements(sob)-1 do begin
           running_snr = [running_snr,max([1.0+10*1./mean(sob[max([0,ii-4]):min([ii+4,n_elements(sob)-1])]/uob[max([0,ii-4]):min([ii+4,n_elements(sob)-1])]),1.5])]
        endfor
;        m = where(wave ge cont_st[i] and wave le cont_en[i] and mob ne 1 and sob gt 0 and sob lt running_snr,mc)
        m = where(wave ge cont_st[i] and wave le cont_en[i] and mob ne 1 and sob gt 0 and sob lt running_snr,mc)
        if mc ne 0 then begin 
            mob [m] = 2
        endif
    endfor

    ;De-select 70% lowest continuum points in each segment using synthesis
    ;Ensure that both ends have continuum points
    for i = 0,nseg-1 do begin
        frac=0.80
        again:
        j = where(wave ge wran[0,i] and wave le wran[1,i],jc)
        temp = smod[j]
        temp = temp[sort(temp)]
        temp = temp[n_elements(temp)*frac]
        k = where(wave ge wran[0,i] and wave le wran[1,i] and smod lt temp and mob eq 2,kc)
        l = where(wave ge wran[0,i] and wave le wran[1,i] and smod ge temp and mob eq 2,lc)
        if  n_elements(where(l lt j[float(jc)*1/3.])) lt 5 or $
            n_elements(where(l gt j[float(jc)*2/3.])) lt 5 then begin
            frac=frac-0.1
            if frac gt 0. then goto,again
        endif
        if kc ne 0 then mob[k] = 0
    endfor

endif

;Avoid strong NLTE-dominated cores in mask
;cmin = min([0.60+(5.0-grav)*0.09,0.97])
;if feh lt -1. then cmin=min([0.65+(5.0-grav)*0.09,0.97])
;if feh lt -2. then cmin=min([0.72+(5.0-grav)*0.09,0.97])
cmin=0.60
if feh lt -1. then cmin=0.65
if feh lt -2. then cmin=0.72

if keyword_set(line_cores) then begin
    for i=0,n_elements(line_cores)-1 do begin
        j=where(abs(wave-line_cores[i]) lt 4. and sob lt cmin,jc)
        if jc ne 0 then mob[j]=0
    endfor
endif

END

PRO make_line,run=run, newsme=newsme

@sme_struct

if keyword_set(newsme) then begin
	line_merge, 'LINELIST/'+line_list, line_atomic, line_species, line_lande, line_depth, line_ref,term_low=line_term_low,term_upp=line_term_upp,short_format=short_format,extra=line_extra,lulande=line_lulande 
endif else begin
	line_merge, 'LINELIST/'+line_list, line_atomic, line_species, line_lande, line_depth, line_ref,term_low=line_term_low,term_upp=line_term_upp,short_format=short_format,extra=line_extra;,lulande=line_lulande 
endelse

nselect=n_elements(line_species)
print,nselect,' spectral lines read from line lists'

if not keyword_set(line_term_low) then line_term_low = strarr(n_elements(line_species))
if not keyword_set(line_term_upp) then line_term_upp = strarr(n_elements(line_species))
if not keyword_set(line_extra   ) then line_extra    = dblarr(3,n_elements(line_species))
if keyword_set(newsme) then begin
	if not keyword_set(line_lulande ) then line_lulande  = dblarr(2,n_elements(line_species))
endif
if not keyword_set(depthmin)      then depthmin      = 0.0d0

; select lines within wavelength segments and deeper than a given depth.
de=0.7
nrline=20.+((8000.-teff)/1.3e3)^4

j=[0]
if run eq 2 then begin
  ;Select all lines in segment
   for i=0,nseg-1 do begin
      k=where(line_atomic[2,*] ge wran[0,i] - de and line_atomic[2,*] le wran[1,i] + de and line_depth ge depthmin,kc)
      if kc ne 0 then j=[j,k]
      ;Always select broad lines when close to them
      if keyword_set(broad_lines) then begin
          for l=0,n_elements(broad_lines)-1 do begin
              if abs(broad_lines[l]-wran[0,i]) lt 100. or abs(broad_lines[l]-wran[1,i]) lt 100. then begin
                  k=where(line_atomic[2,*] eq broad_lines[l],kc)
                  if kc ne 0 then j=[j,k]
              endif
          endfor
      endif
   endfor
endif else begin
   ;Select only strongest lines closest to line mask
    readcol,'DATA/'+line_mask,line0,line_st,line_en,comment=';',/silent,format='d,d,d'
    for i=0,n_elements(line_st)-1 do begin
        k=where(line_atomic[2,*] ge line_st[i] - de and line_atomic[2,*] le line_en[i] + de and line_depth gt depthmin,kc)
        if kc ne 0 then begin 
            k=reverse(k[sort(line_depth[k])])
            j=[j,k[0:min([nrline,kc-1])]]
            ;Always select the main line if it is present
            l=where(line_atomic[2,*] eq line0[i],lc)
            if lc ne 0 then j=[j,l] else print,'Missing in line-list: '+string(line0[i],format='(f10.4)')
            ;Always select broad lines when close to them
            if keyword_set(broad_lines) then begin
                for l=0,n_elements(broad_lines)-1 do begin
                    if abs(broad_lines[l]-line0[i]) lt 100. then begin
                        k=where(line_atomic[2,*] eq broad_lines[l],kc)
                        if kc ne 0 then j=[j,k]
                    endif
                endfor
            endif
        endif
    endfor
endelse

nselect   =  n_elements(j)-1
if nselect gt 0 then j=unique(j[1:*])
nselect   =  n_elements(j)

species   =  line_species  [j]
atomic    =  line_atomic [*,j]
lande     =  line_lande    [j]
depth     =  line_depth    [j]
lineref   =  line_ref      [j]
term_low  =  line_term_low [j]
term_upp  =  line_term_upp [j]
extra     =  line_extra  [*,j]
if keyword_set(newsme) then lulande   =  line_lulande[*,j]

print,nselect,' spectral lines selected within wavelength segments'

;Sort according to wavelength
j     =  sort(atomic[2,*])
species   =  species  [j]
atomic    =  atomic [*,j]
lande     =  lande    [j]
depth     =  depth    [j]
lineref   =  lineref  [j]
term_low  =  term_low [j]
term_upp  =  term_upp [j]
extra     =  extra  [*,j]
if keyword_set(newsme) then lulande   =  lulande[*,j]

END

PRO make_struct,run=run,norm=norm,nlte=nlte,status=status, newsme=newsme

;Generalized routine to create sme input structure. Variables are
;communicated through common blocks in sme_struct.pro

@sme_struct

;Set remaining fixed variables, usually fixed
;version        = 4.37
;if keyword_set(newsme) then version = 5.36
gam6           = 1.0
accwi          = 0.005 
accrt          = 0.005
clim           = 0.01
chirat         = 0.001
nmu            = 7
obs_type       = 3
iptype         = 'gauss'
mu             = reverse(sqrt(0.5*(2*dindgen(nmu)+1)/float(nmu)))
id             = systime()
if glob_free[0] eq '-1' then $
   glob_free      = string(-1)
atmo_pro       = 'interp_atmo_grid' ; 'interpmarcs2012' 
if not keyword_set(atmogrid_file) then $
   atmogrid_file  = 'marcs2012.sav' 
atmo_grid_vers =  4.0
if not keyword_set(ab_free) then $
ab_free        = intarr(99)

; get observations 
make_obs,status
if status eq 0 then goto,error

; set initial RV correction for each segment to zero
if vrad_flag eq 0 then vrad = replicate(0.0,nseg) else vrad=0.

; set initial continuum scale to 1.0 (2 values for each segment, for cscale_flag=1)
if cscale_flag eq 1 then cscale = replicate(1.0,2,nseg)
if cscale_flag eq 0 then cscale = replicate(1.0,nseg)
if cscale_flag lt 0 then cscale = 1.0

; pre-normalise observations if necessary
if keyword_set(norm) and cscale_flag ge 0 then pre_norm

; set the mask
make_mob

; get the line list
make_line,run=run, newsme=newsme

;gf-free
gf_free=lonarr(n_elements(species))
;if keyword_set(lines_gf_free) then begin 
;   for i=0,n_elements(lines_gf_free)-1 do begin
;      select=where(atomic[2,*] eq lines_gf_free[i],nselect)
;      if nselect ne 0 then gf_free[select]=1
;   endfor
;endif

; build new SME structure.


if keyword_set(newsme) then begin
  ;;;;;;;;
  ; The default value for atmo_depth was 'RHOX' until 180212
  ; Changed after email conversation with T. Nordlander to 'TAU'
  ; Changed back to RHOX after email with T. Nordlander with improved interp routines
  ;;;;;;;;
  setdefaultvalue, atmo_depth, 'RHOX';'RHOX'
  setdefaultvalue, atmo_interp, 'TAU'
  setdefaultvalue, atmo_geom, 'PP'
  ; Temporary fix for stagger grid, which has only EITHER tau or rhox depth scale available:
  if strmatch(atmogrid_file, '*stagger-t*') then atmo_depth = 'TAU'
  if strmatch(atmogrid_file, '*stagger-r*') then atmo_interp = 'RHOX'
  
  sme = { $
   version        :  version       , $
   id             :  id            , $
   teff           :  teff          , $
   grav           :  grav          , $
   feh            :  feh           , $
   vmic           :  vmic          , $
   vmac           :  vmac          , $
   vsini          :  vsini         , $
   vrad           :  vrad          , $
   vrad_flag      :  vrad_flag     , $
   cscale         :  cscale        , $
   cscale_flag    :  cscale_flag   , $
   gam6           :  gam6          , $
   accwi          :  accwi         , $
   accrt          :  accrt         , $
   clim           :  clim          , $
   maxiter        :  maxiter       , $
   chirat         :  chirat        , $
   nmu            :  nmu           , $
   abund          :  abund         , $
   mu             :  mu            , $
   atmo : { $
      source :  atmogrid_file , $       
      method :  'grid', $
      depth  :  atmo_depth, $
      interp :  atmo_interp, $
      geom   :  atmo_geom $
   }, $
   sob            :  sob           , $
   uob            :  uob           , $
   obs_name       :  obs_name      , $
   obs_type       :  obs_type      , $
   iptype         :  iptype        , $
   glob_free      :  glob_free     , $
   ab_free        :  ab_free       , $
   gf_free        :  gf_free       , $
   species        :  species       , $
   atomic         :  atomic        , $
   lande          :  lande         , $
   depth          :  depth         , $
   lineref        :  lineref       , $
   line_term_low  :  term_low      , $
   line_term_upp  :  term_upp      , $
   short_format   :  short_format  , $
   line_extra     :  extra         , $
   line_lulande   :  lulande       , $
   nseg           :  nseg          , $
   wran           :  wran          , $
   wave           :  wave          , $
   wind           :  wind          , $
   mob            :  mob           , $
   ipres          :  ipres         , $
   auto_alpha     :  auto_alpha      $
  }
endif else begin
  sme = { $
   version        :  version       , $
   id             :  id            , $
   teff           :  teff          , $
   grav           :  grav          , $
   feh            :  feh           , $
   vmic           :  vmic          , $
   vmac           :  vmac          , $
   vsini          :  vsini         , $
   vrad           :  vrad          , $
   vrad_flag      :  vrad_flag     , $
   cscale         :  cscale        , $
   cscale_flag    :  cscale_flag   , $
   gam6           :  gam6          , $
   accwi          :  accwi         , $
   accrt          :  accrt         , $
   clim           :  clim          , $
   maxiter        :  maxiter       , $
   chirat         :  chirat        , $
   nmu            :  nmu           , $
   abund          :  abund         , $
   mu             :  mu            , $
   atmo_pro       :  atmo_pro      , $
   atmogrid_file  :  atmogrid_file , $       
   atmo_grid_vers :  atmo_grid_vers, $     
   sob            :  sob           , $
   uob            :  uob           , $
   obs_name       :  obs_name      , $
   obs_type       :  obs_type      , $
   iptype         :  iptype        , $
   glob_free      :  glob_free     , $
   ab_free        :  ab_free       , $
   gf_free        :  gf_free       , $
   species        :  species       , $
   atomic         :  atomic        , $
   lande          :  lande         , $
   depth          :  depth         , $
   lineref        :  lineref       , $
   line_term_low  :  term_low      , $
   line_term_upp  :  term_upp      , $
   short_format   :  short_format  , $
   line_extra     :  extra         , $
;   line_lulande   :  lulande       , $
   nseg           :  nseg          , $
   wran           :  wran          , $
   wave           :  wave          , $
   wind           :  wind          , $
   mob            :  mob           , $
   ipres          :  ipres         , $
   auto_alpha     :  auto_alpha      $
  }
endelse

if keyword_set(nlte) then begin
   nltestruct = { $
      nlte_pro          :  'sme_nlte'    , $
      nlte_elem_flags   : nlte_elem_flags, $
      nlte_subgrid_size : [3,3,3,3]      , $
      nlte_grids        : nlte_grids     , $
      nlte_debug        : 1                $
                }
   if keyword_set(newsme) then begin
      sme = create_struct(temporary(sme), 'NLTE', nltestruct)
   endif else begin
      sme = create_struct(temporary(sme), nltestruct)
   endelse
endif


; store sme input structure for reference
save, file='OUTPUT/' + obs_name+'_SME.inp', sme

; check if there are enough points left in line masks before running
i=where(sme.mob eq 1,ic)
if ic lt 5 and (glob_free[0] ne string(-1) or max(ab_free) ne 0) then begin
   print,'Not enough points in line mask, returning'
   status=0
   error:
   return
endif

if keyword_set(run) then begin 
   ;print, 'make_struct: calling sme for star ',obs_name
    sme_main,sme
   ; store sme output structure in idl save file
   save, file='OUTPUT/' + obs_name+'_SME.out', sme

endif

END
