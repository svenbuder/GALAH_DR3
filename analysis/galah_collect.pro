
PRO GALAH_struct,galah

NaN=abs(sqrt(-1))

galah = { $
   FIELD         :  ''        , $
   OBJECT        :  ''        , $
   SOBJECT_ID    :  LONG64(0) , $
   FILENAME      :  ''        , $
   GALAH_ID      :  ''        , $
   FILE_FLAG     :  0.0       , $
   COB_ID        :  ''        , $
   PIVOT         :  ''        , $
   RESOLUTION    :  0.0       , $
   RA            :  0.d0      , $
   DEC           :  0.d0      , $
   EBV           :  0.0       , $
   FLAG_SP       :  0.0       , $
   VEL           :  0.0       , $
   E_VEL         :  0.0       , $
   C_VEL         :  0.0       , $
   TEFF          :  0.0       , $
   E_TEFF        :  0.0       , $
   C_TEFF        :  0.0       , $
   LOGG          :  0.0       , $
   E_LOGG        :  0.0       , $
   C_LOGG        :  0.0       , $
   FE_H          :  0.0       , $
   E_FE_H        :  0.0       , $     
   C_FE_H        :  0.0       , $
   FEH           :  0.0       , $
   E_FEH         :  0.0       , $
   C_FEH         :  0.0       , $
   VMIC          :  0.0       , $
   E_VMIC        :  0.0       , $
   C_VMIC        :  0.0       , $
   VMAC          :  0.0       , $
   E_VMAC        :  0.0       , $
   C_VMAC        :  0.0       , $
   ABFE          :  0.0       , $
   E_ABFE        :  0.0       , $
   ALPHA_FE      :  0.0       , $
   E_ALPHA_FE    :  0.0       , $
   C_ALPHA_FE    :  0.0       , $
   VRAD          :  0.0       , $
   E_VRAD        :  0.0       , $
   C_VRAD        :  0.0       , $
   VSINI         :  0.0       , $
   E_VSINI       :  0.0       , $
   C_VSINI       :  0.0       , $
   TECH          :  ''        , $
   RED_FLAG      :  0.0       , $
   FLAG_GUESS    :  0.0       , $
   FLAG_CANNON   :  0.0       , $
   SNR2_C1_IRAF  :  0.0       , $
   SNR2_C2_IRAF  :  0.0       , $
   SNR2_C3_IRAF  :  0.0       , $
   SNR2_C4_IRAF  :  0.0       , $
   SVEL          :  fltarr(5) , $            
   STEFF         :  fltarr(5) , $
   SFEH          :  fltarr(5) , $
   SLOGG         :  fltarr(5) , $
   SVMIC         :  fltarr(5) , $
   SVSINI        :  fltarr(5) , $
   CHI           :  fltarr(99), $
   ITER          :  lonarr(99), $
   SN            :  fltarr(99), $
   MODE          :  strarr(99), $
   ABUND         :  fltarr(99), $
   A_ABUND       :  fltarr(99), $
   E_ABUND       :  fltarr(99), $
   C_ABUND       :  fltarr(99), $
   AFLAG         :  lonarr(99), $
   ALOG          :  lonarr(99), $
   TIME          :  fltarr(99), $
   LINEFLUX      :  fltarr(99), $
   LINE_FC       :  fltarr(99), $
   BA_VRAD       :  0.0       , $
   LI_VRAD       :  0.0       , $
   SI6722        :  0.0       , $
   E_SI6722      :  0.0       , $
   C_SI6722      :  0.0       , $
   SI_VRAD       :  0.0       , $
   MASS          :  0.0       , $
   AGE           :  0.0       , $
   BC_K          :  0.0       , $
   LBOL          :  0.0         $
}

GALAH.SOBJECT_ID    = NaN
GALAH.FILE_FLAG     = NaN
GALAH.RA            = NaN
GALAH.DEC           = NaN
GALAH.EBV           = NaN
GALAH.RESOLUTION    = NaN
GALAH.RED_FLAG      = NaN
GALAH.FLAG_GUESS    = NaN
GALAH.FLAG_CANNON   = NaN
GALAH.VEL           = NaN
GALAH.E_VEL         = NaN
GALAH.C_VEL         = NaN
GALAH.TEFF          = NaN
GALAH.E_TEFF        = NaN
GALAH.C_TEFF        = NaN
GALAH.LOGG          = NaN
GALAH.E_LOGG        = NaN
GALAH.C_LOGG        = NaN
GALAH.FE_H          = NaN
GALAH.E_FE_H        = NaN
GALAH.C_FE_H        = NaN
GALAH.FEH           = NaN
GALAH.E_FEH         = NaN
GALAH.C_FEH         = NaN
GALAH.VMIC          = NaN
GALAH.E_VMIC        = NaN
GALAH.C_VMIC        = NaN
GALAH.VMAC          = NaN
GALAH.E_VMAC        = NaN
GALAH.ABFE          = NaN
GALAH.E_ABFE        = NaN
GALAH.ALPHA_FE      = NaN
GALAH.E_ALPHA_FE    = NaN
GALAH.C_ALPHA_FE    = NaN
GALAH.VRAD          = NaN
GALAH.E_VRAD        = NaN
GALAH.C_VRAD        = NaN
GALAH.VSINI         = NaN
GALAH.E_VSINI       = NaN
GALAH.C_VSINI       = NaN
GALAH.SVEL[*]       = NaN
GALAH.STEFF[*]      = NaN
GALAH.SFEH[*]       = NaN
GALAH.SLOGG[*]      = NaN
GALAH.SVMIC[*]      = NaN
GALAH.SVSINI[*]     = NaN
GALAH.SNR2_C1_IRAF  = NaN
GALAH.SNR2_C2_IRAF  = NaN
GALAH.SNR2_C3_IRAF  = NaN
GALAH.SNR2_C4_IRAF  = NaN
GALAH.CHI[*]        = NaN
GALAH.SN[*]         = NaN
GALAH.ABUND[*]      = NaN
GALAH.A_ABUND[*]    = NaN
GALAH.E_ABUND[*]    = NaN
GALAH.C_ABUND[*]    = NaN
GALAH.AFLAG[*]      = -1
GALAH.TIME[*]       = NaN
GALAH.LINEFLUX[*]   = NaN
GALAH.LINE_FC[*]    = NaN
GALAH.BA_VRAD       = NaN
GALAH.LI_VRAD       = NaN
GALAH.SI6722        = NaN
GALAH.E_SI6722      = NaN
GALAH.C_SI6722      = NaN
GALAH.SI_VRAD       = NaN
GALAH.AGE[*]        = NaN
GALAH.MASS[*]       = NaN
GALAH.BC_K[*]       = NaN
GALAH.LBOL[*]       = NaN

END

PRO GALAH_collect,final,dir=dir,sp=sp,to=to,lo=lo,fo=fo,offset_free=offset_free,offset_lbol=offset_lbol,silent=silent

;print,'Running offset_lbol true and silent true'
;offset_lbol = 1
;silent = 1

  print,'You have the options sp, to, lo, fo, offset_free, offset_lbol with'
  print,'sp == only collect parameters'
  print,'to == Teff offset'
  print,'lo == logg offset'
  print,'fo == feh offset'
  print,'offset_free == to=0,lo=0.15,fo=0.1'
  if keyword_set(offset_free) then begin
     lo = 0.15
     fo = 0.1
  endif
  print,'offset_lbol == to=0,lo=0,fo=0.1' 
  if keyword_set(offset_lbol) then fo=0.1

;print,'fo',fo
  
;ON_IOERROR,next
cla=command_line_args(count=count)
if count ne 0 then begin
   final  = cla[0]
endif

get_date,dte
if not keyword_set(dir) then dir='OUTPUT'

;Read list of object and elements
readcol,final,field,object,setup,format='a,a,a',/silent,comment=';'

openw,good_sp,final+'_NoTech',/get_lun

;Read iraf_reduction input fits
input = mrdfits('DATA/sobject_iraf_53_2MASS_GaiaDR2_WISE_PanSTARRSDR1_BailerJones_K2seis_small.fits',1,/silent)
cannon = mrdfits('DATA/sobject_iraf_iDR2_cannon_small.fits',1,/silent)

;Solar abundances
marcs_abund,sun,elem,fehmod=0.0,eonh12=sunH12

;Read template fits-files
GALAH_struct,res

res = replicate(res,n_elements(object))

;Chisq-SN relation
x  = [ 0. ,1  ,1.5,1.8,2.0,2.4,2.6]
y  = [-1.1,1.2,1.4,1.8,2.2,3.6,5.0]
xchi = lindgen(101)*0.03
ychi = interpol(y,x,xchi,/spline)

for i=0L,n_elements(object)-1 do begin

   res[i].object     = object[i]
   res[i].sobject_id = long64(object[i])
   res[i].field      = field[i]
   if not keyword_set(silent) then print,object[i],res[i].sobject_id

   if i mod 50L eq 0 then  print,float(i)/float(n_elements(object))*100,'%  ',systime(),format='(f10.3,a3,a-35)'

   ;fits='SPECTRA/irafdr51/'+strmid(object[i],0,6)+'/combined/'+object[i]+'1.fits'
   fits='SPECTRA/dr5.3/'+strmid(object[i],0,6)+'/standard/com/'+object[i]+'1.fits'
   fits_exists = file_test(fits)
   if not fits_exists then begin
      fits='SPECTRA/dr5.3/'+strmid(object[i],0,6)+'/standard/com2/'+object[i]+'1.fits'
      fits_exists = file_test(fits)
      if not fits_exists then begin
         fits='SPECTRA/dr5.2/'+strmid(object[i],0,6)+'/standard/com/'+object[i]+'1.fits'
         fits_exists = file_test(fits)
      endif
   endif

   if fits_exists then begin

      ext0 = mrdfits(fits,0,h0,/silent)
      fits_info,fits,/silent,n_ext=n_ext
      res[i].filename = str_replace(fits,'SPECTRA/'+field[i]+'/','') 
      res[i].resolution = 0.0
      high_low = fxpar(h0,'SLITMASK')
      if high_low eq 'IN      ' then res[i].resolution = 1.0

      ;find object in iraf reduction input
      j=where(object[i] eq input.sobject_id,jc) & j=j[0]
      
      if jc eq 0 then begin
         print,'object not in sobject_iraf_53 file'
         return
      endif

      ;res[i].galah_id     = input[j].galah_id
      res[i].ra           = input[j].ra
      res[i].dec          = input[j].dec
      ;print,res[i].ra,res[i].dec
      res[i].ebv          = input[j].ebv
      ;res[i].snr2_c1_iraf = input[j].snr_c1_iraf
      res[i].snr2_c2_iraf = input[j].snr_c2_iraf
      ;res[i].snr2_c3_iraf = input[j].snr_c3_iraf
      ;res[i].snr2_c4_iraf = input[j].snr_c4_iraf
      res[i].red_flag     = input[j].red_flag
      res[i].flag_guess   = input[j].flag_guess

      jj=where(object[i] eq cannon.sobject_id,jjc) & jj=jj[0]

      if jjc gt 0 then begin
         res[i].flag_cannon  = cannon[jj].flag_cannon
      endif else begin
         res[i].flag_cannon  = -1
      endelse

   endif else if not keyword_set(silent) then $
      print,'No FITS found for ',object[i]


   up2:
   readcol,'mode_'+setup[i]+'_collect',modeall,l0,ms,me,format='a,d,d,d',/silent,comment=';'
   mode=unique(modeall)
   res.mode[0:n_elements(mode)-1]=mode
   if keyword_set(sp) then begin 
      if sp eq 'Sp' then mode='Sp' $ ;Collect only stellar parameters from Sp-files
      else mode[0]=sp
   endif
   
   for ii=0,n_elements(mode)-1 do begin

      file = 'SMALL_OUTPUT/'+mode[ii]+'/OUTPUT_'+field[i]+'/'+field[i]+'_'+object[i]+'_'+setup[i]+'_'+mode[ii]+'_SME.out'
      log  = 'SMALL_OUTPUT/'+mode[ii]+'/OUTPUT_'+field[i]+'/idl_'+field[i]+'_'+object[i]+'_'+setup[i]+'_'+mode[ii]+'.log'

      file_exists = file_test(file)
      log_exists = file_test(log)

      by_karin = ['Li6708','Li6708_NoRV','C6588','O','Na','Mg5711','Al','Si','K7699','Ca','Sc','V4832','Cr','Mn','Ba']
      was_run_by_karin = where(mode[ii] eq by_karin,run_by_karin)
      if run_by_karin eq 0 then begin
         if not file_exists then file=dir+'/'+field[i]+'_'+object[i]+'_'+setup[i]+'_'+mode[ii]+'_SME.out'
         if not log_exists then log=dir+'/idl_'+field[i]+'_'+object[i]+'_'+setup[i]+'_'+mode[ii]+'.log'
      endif

      if mode[ii] eq 'Li6708_NoRV' then begin
         file = 'SMALL_OUTPUT/Li6708_NoRV/'+field[i]+'_'+object[i]+'_'+setup[i]+'_Li6708_SME.out'
         log  = 'SMALL_OUTPUT/Li6708_NoRV/idl_'+field[i]+'_'+object[i]+'_'+setup[i]+'_Li6708.log'
         mode[ii] = 'Li6708' ; to actually collect the values
      endif

      file_exists = file_test(file)
      log_exists = file_test(log)

      if log_exists then if file_lines(log) ne 0 then begin 

         res[i].alog[ii]=1

         openr,unit,log,/get_lun
         logtext=strarr(file_lines(log))
         readf,unit,logtext
         free_lun,unit

         if ii eq 0 then begin
            timeiter=where(strpos(logtext,'Finishing') eq 0,i0) & timeiter=timeiter[-1]
            if timeiter ne -1 then time=float(strsplit(strmid(logtext[timeiter],24),'hms',/extract))
         endif else begin
            iter=where(strpos(logtext,'Finally') eq 0,i0) & iter=iter[-1]
            if iter ne -1 then time=float(strsplit(strmid(logtext[iter],33),'hms',/extract))
         endelse
         if i0 ne 0 then begin 
            if i0 gt 1 then begin
               iter=iter[-1]
               if not keyword_set(silent) then print,'double entry for '+string(object[i])
            endif
            res[i].time[ii]=time[0]*60*60.+time[1]*60.+time[2]
         endif

         if ii eq 0 then begin 
            iter=where(strpos(logtext,'Input  '+mode[0]) ne -1,i0) & iter=iter[-1]
            if i0 ne 0 then begin 
               dum = strsplit(logtext[iter],' ',/extract)
               res[i].steff[0] =  dum[1]
               res[i].slogg[0] =  dum[2]
               res[i].sfeh[0]  =  dum[3]
               res[i].svel[0]  =  dum[4]
            endif
            iter=where(strpos(logtext,'Initial  '+mode[0]) ne -1,i1)
            if i1 ne 0 then begin 
               dum = strsplit(logtext[iter],' ',/extract)
               res[i].steff[1] =  dum[1]
               res[i].slogg[1] =  dum[2]
               res[i].sfeh[1]  =  dum[3]
               res[i].svel[1]  =  dum[4]
               res[i].vel      =  dum[4]
            endif
            iter=where(strpos(logtext,'2  '+mode[0]) ne -1,i2)
            if i2 ne 0 then begin 
               dum = strsplit(logtext[iter],' ',/extract)
               res[i].steff[2] =  dum[1]
               res[i].slogg[2] =  dum[2]
               res[i].sfeh[2]  =  dum[3]
               res[i].svmic[2] =  dum[5]
               res[i].svsini[2]=  dum[6]
               res[i].svel[2]  =  dum[7]
               res[i].vel      =  res[i].vel + dum[7]
            endif
            iter=where(strpos(logtext,'4  '+mode[0]) ne -1,i4)
            if i4 ne 0 then begin 
               dum = strsplit(logtext[iter],' ',/extract)
               res[i].steff[3] =  dum[1]
               res[i].slogg[3] =  dum[2]
               res[i].sfeh[3]  =  dum[3]
               res[i].svmic[3] =  dum[5]
               res[i].svsini[3]=  dum[6]
               res[i].svel[3]  =  dum[7]
               res[i].vel      =  res[i].vel + dum[7]
            endif
            iter=where(strpos(logtext,'6  '+mode[0]) ne -1,i6)
            if i6 ne 0 then begin 
               dum = strsplit(logtext[iter],' ',/extract)
               res[i].steff[4] =  dum[1]
               res[i].slogg[4] =  dum[2]
               res[i].sfeh[4]  =  dum[3]
               res[i].svmic[4] =  dum[5]
               res[i].svsini[4]=  dum[6]
               res[i].svel[4]  =  dum[7]
               res[i].vel      =  res[i].vel + dum[7]
            endif
            iter=where(strpos(logtext,'Interim values: GRAV_OUT') ne -1,isp)
            if isp ne 0 then begin
               dum = strsplit(logtext[iter[-1]],' ',/extract)
               res[i].mass     =  strsplit(dum[5],',',/extract)
               res[i].age      =  strsplit(dum[7],',',/extract)
               res[i].lbol     =  dum[9]
            endif
            iter=where(strpos(logtext,'BC_K') ne -1, ibc)
            if ibc ne 0 then begin
               dum = strsplit(logtext[iter[-1]],' ',/extract)
               res[i].bc_k     = dum[6]
            endif

            ; Checking for gridlimit
            iter=where(strpos(logtext,'returning.') ne -1,igrid)
            if igrid ne 0 then res[i].flag_sp += 2

            ; Checking for Gaussian RV fail
            iter=where(strpos(logtext,'Gaussian fit at Balmer lines') ne -1,igauss)
            if igauss ne 0 then res[i].flag_sp += 4

            ; Checking for fail of ELLI age estimation
            iter=where(strpos(logtext,'no sys.argv') ne -1,ielli)
            if ielli ne 0 then res[i].flag_sp += 8

            ; Check for time-out on ISAAC
            if res[i].flag_sp eq 0 and timeiter eq -1 then res[i].flag_sp += 16

            if not keyword_set(silent) then print,'flag_sp: ',res[i].flag_sp

         endif
      endif else begin
         print,'No LOG found for ',object[i]
      endelse

      if file_exists then begin

         ; Because sme_clean_smod
         sme_clean_smod = [-99]

         restore,file

         if (sme.glob_free(0) ne '      -1' or max(sme.ab_free) eq 1) and $
            sme.maxiter eq 20 then begin
               inf  = where(finite([sme.teff,sme.grav,sme.feh,sme.vsini,sme.abund,sme.vrad]) eq 0,infc)
               ;for u=0,n_elements(sme.pname)-1 do sme.punc[u]=sqrt(sme.covar[u,u])
            endif else infc = 1

         ; compute mean squared difference of full and clean smod within final
         ; line mask 
         if sme_clean_smod[0] ne -99 then begin
            final_mob = where(sme.mob eq 1)
            res[i].line_fc[ii] = mean((sme_full_smod[final_mob]-sme_clean_smod[final_mob])^2)
         endif

         if ii eq 0 then begin

            if infc then res[i].tech='conv|' ;Not converged
            if infc then res[i].flag_sp += 1

            if infc eq 0 then begin 
               line_mask  = where(sme.mob eq 1)
               res[i].iter[ii] = n_elements(sme.rchisq)
               res[i].chi[ii]  = sme.rchisq(res[i].iter[ii]-1)  ;total((sme.smod-sme.sob)^2/sme.uob^2)/(n_elements(sme.wave)-2-1.) ;
               res[i].sn[ii]   = mean(sme.sob[line_mask]/abs(sme.uob[line_mask]))

               jt = where(sme.pname eq 'TEFF')              & jt=jt[0]
               jg = where(sme.pname eq 'GRAV')              & jg=jg[0]
               jf = where(sme.pname eq 'FEH')               & jf=jf[0]
               js = where(sme.pname eq 'VSINI')             & js=js[0]
               jr = where(sme.pname eq 'VRAD')              & jr=jr[0]
               jc = where(strpos(sme.pname,'ABUND') ne -1 ) & jc=jc[0]
            endif else begin 
               jt=-1 & jg=-1 & jf=-1 & js=-1 & jr=-1 & jc=-1
            endelse

            res[i].teff    = sme.teff
            if keyword_set(to) then res[i].teff = sme.teff+to
            if jt ne -1 then begin 
               res[i].e_teff  = sme.punc[jt]
               res[i].c_teff  = sqrt(sme.covar[jt,jt])
            endif
            if not keyword_set(silent) then print,'teff', res[i].teff,res[i].e_teff,res[i].c_teff
            res[i].logg    = sme.grav
            if keyword_set(lo) then res[i].logg = sme.grav+lo
            if strmid(field[i],strlen(field[i])-10,10) eq 'logg_fixed' then begin
               res[i].logg = alog10(numax[i]/3090.0*sqrt(res[i].teff/5777.0)*10.0^4.43770565)
               ;print,object[i],' numax: ',numax[i],' and teff ',res[i].teff,' -> seismic logg ',res[i].logg  
            endif
            if jg ne -1 then begin 
               res[i].e_logg  = sme.punc[jg]
               res[i].c_logg  = sqrt(sme.covar[jg,jg])
            endif
            res[i].feh     = sme.feh
            if keyword_set(fo) then res[i].feh = sme.feh+fo
            if jf ne -1 then begin 
               res[i].e_feh  = sme.punc[jf]
               res[i].c_feh  = sqrt(sme.covar[jf,jf])
            endif
            if not keyword_set(silent) then print,'feh', res[i].feh,res[i].e_feh,res[i].c_feh
            res[i].vsini   = sme.vsini
            if js ne -1 then begin 
               res[i].e_vsini = sme.punc[js]
               res[i].c_vsini = sqrt(sme.covar[js,js])
               if finite(sme.punc[js]) and finite(res[i].teff) then printf,good_sp,field[i],object[i],setup[i],format='(a-20,a-30,a-10)'
            endif
            res[i].vrad    = sme.vrad
            res[i].vel     = res[i].vel + sme.vrad
            if jr ne -1 then begin 
               res[i].e_vrad  = sme.punc[jr]
               res[i].e_vel   = sme.punc[jr]
               res[i].c_vrad  = sqrt(sme.covar[jr,jr])
               res[i].c_vel   = sqrt(sme.covar[jr,jr])
            endif
            if jc ne -1 then begin 
               dum            = where(sme.ab_free) & dum=dum[0]
               res[i].abfe    = sme.abund[dum]-sun[dum]
            endif

            res[i].vmic = sme.vmic
            res[i].vmac = sme.vmac

            ;if finite(sme.punc[js]) then printf,good_sp,field[i],object[i],setup[i],format='(a-20,a-30,a-10)'
               
         endif else begin
            
            if max(sme.ab_free) eq 1 then begin 

               ml=strlen(mode[ii])
               if ml gt 2 then ml=ml-4

               if mode[ii] eq 'Ba5854' then res[i].ba_vrad = sme.vrad
               if mode[ii] eq 'Li6708' and min(where(fs(sme.glob_free) eq 'VRAD')) ne -1 then res[i].li_vrad = sme.vrad
               if mode[ii] eq 'Si6722' then begin
                  res[i].si_vrad = sme.vrad
                  res[i].si6722 = sme.abund[where(sme.ab_free)]-sun[where(sme.ab_free)]
                  res[i].e_si6722 = sme.punc[-1]
                  res[i].c_si6722 = sqrt(sme.covar[-1,-1])
               endif

               correct_structure_format = tag_exist(sme, 'RCHISQ')
               if not correct_structure_format then goto,next

               line_mask  = where(sme.mob eq 1)
               res[i].iter[ii] = n_elements(sme.rchisq)
               res[i].chi[ii]  = sme.rchisq(res[i].iter[ii]-1)
               res[i].sn[ii]   = mean(sme.sob[line_mask]/abs(sme.uob[line_mask]))

               ;Check significance of line-depth segment by segment
               detect=fltarr(sme.nseg)-1
               for iii=0,sme.nseg-1 do begin 
                  seg=where(sme.wave ge sme.wran[0,iii] and sme.wave le sme.wran[1,iii])
;                  cont=mean([sme.smod[seg[0]],sme.smod[seg[n_elements(seg)-1]]]) ;Estimate of model continuum in segment
                  cont=mean(sme_full_smod[where(sme_full_mob eq 2)])
                  line_mask_seg  = where(sme.mob eq 1 and sme.wave ge sme.wran[0,iii] and sme.wave le sme.wran[1,iii])
                  if line_mask_seg[0] ne -1 then begin 
                     sn   = mean(sme.sob[line_mask_seg]/abs(sme.uob[line_mask_seg]))
                     ;if cont-min(sme.smod[line_mask_seg]) gt max([1./sn,0.02]) then detect[iii]=0 ;At least upper limit
                     if cont-min(sme.smod[line_mask_seg]) gt max([1.5/sn,0.03]) then detect[iii]=1 ;Detection
                     if not finite(res[i].lineflux[ii]) then res[i].lineflux[ii] = cont-min(sme.smod[line_mask_seg])
                     if cont-min(sme.smod[line_mask_seg]) gt res[i].lineflux[ii] then res[i].lineflux[ii] = cont-min(sme.smod[line_mask_seg])
;                     res[i].lineflux[ii] =
;                     max([res[i].lineflux[ii],cont-min(sme.smod[line_mask_seg])])
                     if cont-min(sme.smod[line_mask_seg]) lt max([1.5/sn,0.03]) then res[i].lineflux[ii]= max([1.5/sn,0.03])
                  endif
               endfor
               ;print,mode[ii],detect
               ;print,mode[ii]
               ;print,res[i].sn[ii]


               ;if mode[ii] eq 'Ba5854' then begin
               ;   print,mode[ii]
               ;   print,max(sme.smod)-min(sme.smod[line_mask]),' gt',max([1./res[i].sn[ii],0.05])
               ;   print,total(sme.depth[where(strmid(sme.species,0,ml) eq strmid(mode[ii],0,ml))]),' gt ',0.3
               ;   print,finite(sme.punc[0]),' True==1'
               ;   print,sme.punc[-1],' lt ',0.3
               ;   print,n_elements(line_mask),' ge ',3
               ;   print,max(detect)
               ;endif

               if mode[ii] eq 'Fe' then begin
                  if not keyword_set(silent) then print,'fe'
                  if not keyword_set(silent) then print,sme.abund[25]
                  res[i].fe_h   = sme.abund[25] - 10^(sme.abund[0]) - (sun[25] - 10^(sun[0])) + res[i].feh
                  res[i].e_fe_h = sqrt((sme.punc[-1])^2. + (res[i].e_feh)^2.)
                  res[i].c_fe_h = sqrt((sqrt(sme.covar[-1,-1]))^2. + (res[i].e_feh)^2.)
               endif
               
               if mode[ii] eq 'Li6708' and min(where(fs(sme.glob_free) eq 'VRAD')) ne -1 then begin 
                  minli=where(sme.smod eq min(sme.smod) and sme.mob eq 1) ;if minimum depth is outside line mask, reject
                  if minli[0] eq -1 then begin 
                     detect=0
                     print,'Li6708 rejected ',res[i].field,'     ',res[i].sobject_id
                  endif
               endif

               ;print,sme.punc[-1]
               ;print,n_elements(line_mask)
               ;print,mode[ii],max(detect),sme.abund[where(sme.ab_free)]-sun[where(sme.ab_free)]
               
               if max(detect) ge 0 and $
                  total(sme.depth[where(strmid(sme.species,0,ml) eq strmid(mode[ii],0,ml))]) gt 0.3 and $
                  finite(sme.abund[where(sme.ab_free)]) and $
                  finite(sme.punc[-1]) and $
                  ;sme.punc[-1] lt 0.3 and $
                  n_elements(line_mask) ge 3. then begin

                  if max(detect) eq 1 then begin
                     res[i].abund[ii]   = sme.abund[where(sme.ab_free)]-sun[where(sme.ab_free)]
                     res[i].a_abund[ii] = (sme_abundances(sme))[where(sme.ab_free)]
                     res[i].e_abund[ii] = sme.punc[-1]
                     res[i].c_abund[ii] = sqrt(sme.covar[-1,-1])
                     res[i].aflag[ii]   = 0
                  endif else begin 
                     res[i].abund[ii]  = sme.abund[where(sme.ab_free)]-sun[where(sme.ab_free)]
                     res[i].a_abund[ii] = (sme_abundances(sme))[where(sme.ab_free)]
                     res[i].e_abund[ii] = sme.punc[-1]
                     res[i].c_abund[ii] = sqrt(sme.covar[-1,-1])
                     res[i].aflag[ii]  = 1
                  endelse
                  if alog10(res[i].chi[ii])-interpol(ychi,xchi,alog10(res[i].sn[ii])) gt 0 then res[i].aflag[ii]=2
               endif
               if not keyword_set(silent) then print,mode[ii],max(detect),res[i].abund[ii],res[i].e_abund[ii],res[i].c_abund[ii],res[i].aflag[ii]
            endif
         endelse
      endif else begin
         if not keyword_set(silent) then print,'No SME OUT found for ',object[i],' in mode ',mode[ii]
      endelse
      next:
   endfor

endfor

ii=where(mode eq 'Mg' or mode eq 'Si' or mode eq 'Ti')
for i=0,n_elements(object)-1 do begin
   j=where(finite(res[i].abund[ii]) and res[i].aflag[ii] eq 0,jc)
   if jc ne 0 then begin
      res[i].alpha_fe   = total(res[i].abund[ii[j]]/res[i].e_abund[ii[j]]^2)/total(1./res[i].e_abund[ii[j]]^2)
      res[i].e_alpha_fe = sqrt(1./total(1./res[i].e_abund[ii[j]]^2))
      ;print,res[i].alpha_fe,res[i].e_alpha_fe,res[i].e_teff
   endif
endfor

mwrfits,res,final+'.fits',/create

free_lun,good_sp

END

