PRO galah_limits,sub,debug=debug

subset='GALAH_'+fs(sub)+'_lbol'
print,subset

if keyword_set(debug) then initial = mrdfits(subset+'.fits',1)
t = mrdfits(subset+'.fits',1)
u = mrdfits(subset+'.fits',1)

vsini = t.vsini
bad_vsini = where(t.vsini gt 199.)
vsini[bad_vsini] = 199.

readcol,'mode_DR3_collect',modeall,l0,ms,me,format='a,d,d,d',/silent,comment=';'
mode=unique(modeall)

;elstr=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','']
elstr = [ 'H', 'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne' $
       ,'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca' $
       ,'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn' $
       ,'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr' $
       ,'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn' $
       ,'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd' $
       ,'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb' $
       ,'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg' $
       ,'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th' $
       ,'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cs', 'Es' ]

run_upper = $
['Li6708','C6588','O','Na','Mg5711','Al','Si','K7699','Ca','Sc',$
 ;not run because not useful 'Ti4782','Ti4798','Ti4802','Ti5739','Ti5866','Ti6599','Ti6717','Ti7853','Ti4720','Ti4765','Ti4799','Ti4849','Ti4866',
 ;not run because not useful 'Co5647','Co6490','Co6551',;'Mo5751','Mo5858','Sm4837','Sm4854', $
 'Ti4758','Ti4759','Ti4778','Ti4820','Ti5689','Ti5716','Ti5720','Ti4874',$
 'V4784','V4797','V4832','Cr','Mn','Co4781','Co4900','Co6632','Co6679','Co7713','Co7838',$
 'Ni5847','Ni6586','Cu5700','Cu5782','Zn4722','Zn4811','Rb7800',$
 'Sr6550','Y','Y4820','Y4855','Y4884','Y5663','Y5729','Zr4739','Zr4772','Zr4806','Zr4828','Zr5681',$
 'Mo5689','Mo6619','Ru4758','Ru4869','Ru5699',$
 'Ba','La4716','La4749','La4804','La5806','Ce4774',$
 'Nd4811','Nd5741','Nd5770','Nd5812','Nd5842','Sm4720','Sm4792','Sm4848','Eu5819','Eu6645']

for mode_index=0,n_elements(mode)-1 do begin

   if strlen(mode[mode_index]) gt 2 then elem=strmid(mode[mode_index],0,strlen(mode[mode_index])-4) else elem=mode[mode_index]
   atom=where(elstr eq elem)   ;Chose element from list                                                                                                                                                                  
   atom=atom[0]

   print,elem,atom

   do_run = where(mode[mode_index] eq run_upper)

   in_t = where(fs(t[0].mode) eq mode[mode_index]) & in_t=in_t[0]

   if in_t eq -1 then begin
      print,'Element '+mode[mode_index]+' not in t[0].mode!'
      goto,end_element
   endif

   if mode[mode_index] eq 'Li6708' then begin
      print,'Merging Li6708 and Li6708_NoRV'
      
      ; find spectra, where we want to use the Li6708_NoRV measurements
      ; This is the opposite of Karin's approach to find good Li6708_RV ones
      i=where(t.aflag[in_t] ne 0 or t.lineflux[in_t] le 0.15)

      no_rv = mrdfits('SMALL_OUTPUT/Li6708_NoRV/'+subset+'_NLTE.fits',1)

      t[i].CHI[in_t] = no_rv[i].CHI[1]
      t[i].ITER[in_t] = no_rv[i].ITER[1]
      t[i].SN[in_t] = no_rv[i].SN[1]
      t[i].ABUND[in_t] = no_rv[i].ABUND[1]
      t[i].A_ABUND[in_t] = no_rv[i].A_ABUND[1]
      t[i].E_ABUND[in_t] = no_rv[i].E_ABUND[1]
      t[i].C_ABUND[in_t] = no_rv[i].C_ABUND[1]
      t[i].AFLAG[in_t] = no_rv[i].AFLAG[1]
      t[i].ALOG[in_t] = no_rv[i].ALOG[1]
      t[i].TIME[in_t] = no_rv[i].TIME[1]
      t[i].LINEFLUX[in_t] = no_rv[i].LINEFLUX[1]
      t[i].LINE_FC[in_t] = no_rv[i].LINE_FC[1]

      u[i].CHI[in_t] = no_rv[i].CHI[1]
      u[i].ITER[in_t] = no_rv[i].ITER[1]
      u[i].SN[in_t] = no_rv[i].SN[1]
      u[i].ABUND[in_t] = no_rv[i].ABUND[1]
      u[i].A_ABUND[in_t] = no_rv[i].A_ABUND[1]
      u[i].E_ABUND[in_t] = no_rv[i].E_ABUND[1]
      u[i].C_ABUND[in_t] = no_rv[i].C_ABUND[1]
      u[i].AFLAG[in_t] = no_rv[i].AFLAG[1]
      u[i].ALOG[in_t] = no_rv[i].ALOG[1]
      u[i].TIME[in_t] = no_rv[i].TIME[1]
      u[i].LINEFLUX[in_t] = no_rv[i].LINEFLUX[1]
      u[i].LINE_FC[in_t] = no_rv[i].LINE_FC[1]

      ;t[i].CHI[in_t] = t[i].CHI[in_t+1]
      ;t[i].ITER[in_t] = t[i].ITER[in_t+1]
      ;t[i].SN[in_t] = t[i].SN[in_t+1]
      ;t[i].ABUND[in_t] = t[i].ABUND[in_t+1]
      ;t[i].A_ABUND[in_t] = t[i].A_ABUND[in_t+1]
      ;t[i].E_ABUND[in_t] = t[i].E_ABUND[in_t+1]
      ;t[i].C_ABUND[in_t] = t[i].C_ABUND[in_t+1]
      ;t[i].AFLAG[in_t] = t[i].AFLAG[in_t+1]
      ;t[i].ALOG[in_t] = t[i].ALOG[in_t+1]
      ;t[i].TIME[in_t] = t[i].TIME[in_t+1]
      ;t[i].LINEFLUX[in_t] = t[i].LINEFLUX[in_t+1]
      ;t[i].LINE_FC[in_t] = t[i].LINE_FC[in_t+1]

      ;u[i].CHI[in_t] = t[i].CHI[in_t+1]
      ;u[i].ITER[in_t] = t[i].ITER[in_t+1]
      ;u[i].SN[in_t] = t[i].SN[in_t+1]
      ;u[i].ABUND[in_t] = t[i].ABUND[in_t+1]
      ;u[i].A_ABUND[in_t] = t[i].A_ABUND[in_t+1]
      ;u[i].E_ABUND[in_t] = t[i].E_ABUND[in_t+1]
      ;u[i].C_ABUND[in_t] = t[i].C_ABUND[in_t+1]
      ;u[i].AFLAG[in_t] = t[i].AFLAG[in_t+1]
      ;u[i].ALOG[in_t] = t[i].ALOG[in_t+1]
      ;u[i].TIME[in_t] = t[i].TIME[in_t+1]
      ;u[i].LINEFLUX[in_t] = t[i].LINEFLUX[in_t+1]
      ;u[i].LINE_FC[in_t] = t[i].LINE_FC[in_t+1]

   endif

   if do_run ne -1 then begin

      ;the upper grids have:
      ;lineflux_all,ag,tg,gg,fg,vg
      restore,'upper/upper_'+mode[mode_index]+'.sav'

      for i=0,n_elements(t)-1 do begin
         
         if t[i].teff gt 3800 and t[i].flag_sp eq 0 and t[i].feh lt 1 and t[i].lineflux[in_t] lt 1 then begin 

            ti=where(tg lt t[i].teff) & ti=[max(ti),max(ti)+1] ;& print,'Teff  :',tg[ti]
            gi=where(gg lt t[i].logg) & gi=[max(gi),max(gi)+1] ;& print,'logg  :',gg[gi]
            fi=where(fg lt t[i].feh)  & fi=[max(fi),max(fi)+1] ;& print,'[Fe/H]:',fg[fi]   
            vi=where(vg lt vsini[i])& vi=[max(vi),max(vi)+1] ;& print,'Vsini :',vg[vi]   
            lflux=lineflux_all[*,ti,gi,fi,vi]
            tflux=lflux[*,0,*,*,*]+(lflux[*,1,*,*,*]-lflux[*,0,*,*,*])*(t[i].teff -tg[ti[0]])/(tg[ti[1]]-tg[ti[0]])
            gflux=tflux[*,0,0,*,*]+(tflux[*,0,1,*,*]-tflux[*,0,0,*,*])*(t[i].logg -gg[gi[0]])/(gg[gi[1]]-gg[gi[0]])
            fflux=gflux[*,0,0,0,*]+(gflux[*,0,0,1,*]-gflux[*,0,0,0,*])*(t[i].feh  -fg[fi[0]])/(fg[fi[1]]-fg[fi[0]])
            vflux=fflux[*,0,0,0,0]+(fflux[*,0,0,0,1]-fflux[*,0,0,0,0])*(vsini[i]-vg[vi[0]])/(vg[vi[1]]-vg[vi[0]])
            ;plot,ag,vflux,/ylog
            delt=vflux[1:*]-vflux[0:n_elements(ag)-2]
            ii=where(delt gt 0)
            ab=interpol(ag[ii],alog10(vflux[ii]),alog10(t[i].lineflux[in_t])) 

            marcs_abund,abund,fehmod=t[i].feh
            abund[2:*] += t[i].feh
            abund[1:*] = 10^abund[1:*]
            a_x = alog10(abund / abund[0]) + 12
            ab_x = interpol(a_x[atom] + ag[ii],alog10(vflux[ii]),alog10(t[i].lineflux[in_t]))

            ;oplot,ag[ii],vflux[ii],col=2
            ;oplot,[ab],[r[i].lineflux],psym=8,col=4,symsize=3
            ;oplot,[r[i].abund],[r[i].lineflux],psym=8,col=3,symsize=3
            ;read,ok
            if ab lt 4. then begin 
               u[i].abund[in_t]=ab
               u[i].a_abund[in_t]=ab_x
               ;If there is no measurement, adopt upper limit
               ;if upper limit is much more conservative and
               ;lineflux is low then adopt it also
               if (t[i].aflag[in_t] eq -1) or (t[i].aflag[in_t] eq 0 and ab-t[i].abund[in_t] gt 0.2 and t[i].lineflux[in_t] lt 0.08) then begin 
                  t[i].abund[in_t]=ab
                  t[i].a_abund[in_t]=ab_x
                  t[i].aflag[in_t]=1
               endif
            endif

            if keyword_set(debug) then begin
               print,'[X/Fe]'
               print,[t[i].feh,initial[i].abund[in_t],initial[i].aflag[in_t],t[i].abund[in_t],t[i].aflag[in_t],u[i].abund[in_t]]
               print,'A(X)'
               print,[t[i].feh,initial[i].a_abund[in_t],initial[i].aflag[in_t],t[i].a_abund[in_t],t[i].aflag[in_t],u[i].a_abund[in_t]]
            endif
         endif

      endfor
      
   endif else begin
      print,'Will not run upper limits for '+mode[mode_index]                                                                                                                                                      
   endelse

end_element:

endfor

mwrfits,u,subset+'_upper.fits',/create
mwrfits,t,subset+'_final.fits',/create

END
