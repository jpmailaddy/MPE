;calculateBLC.pro
;created by Ethan D. Peck
;This program will contain procedures and functions to calculate the best fit sine curve to the 0 and 90 degree MEPED data.
;Then it will take the intersection point of the BLC angle and the sine curve to get a flux. That flux will be assumed for all angles of the BLC.
;The BLC is calculated by taking the "integral" of the BLC boundary flux and multiplying it by the BLC boundary angle (in radians).

function advanceDate,cdoy
  ;useful function that advances date
  cdoy = cdoy + 1.
  CALDAT,cdoy,dmo,ddy,dyr
  date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]
  return, date
end


pro fsin_old, X, A, F, pder
  ;calculates sine function with 2 parameters for the form y = Asin(Bx)
  ;X is the independent variable
  ;A is an array of free paramenters such that: A[0] = A and A[1] = B
  ;F is the function return, thus F = y.
  ;pder is the derivative of F with respect to each free parameter in A[].
  F = A[0]*sin(A[1]*X)
  pder = [sin(A[1]*X),A[0]*X*cos(A[1]*X)]
end

pro fsin, X, A, F, pder
  ;calculates sine function with 1 parameters for the form y = Asin(x)
  ;X is the independent variable
  ; X = findgen(91)
  ; X = X*!pi/180.
  ;A is an array of free paramenters such that: A[0] = A and A[1] = B
  ;F is the function return, thus F = y.
  ;pder is the derivative of F with respect to each free parameter in A[].
  ;print, 'A = ', A
  ;  stop
  F = A[0]*sin(X)
  pder = [sin(X)]
  ;    print, 'F = ', F
  ;    print, 'dF/dA = ', pder
   ;  stop
end

pro fcos, X, A, F, pder
  ;calculates sine function with 1 parameters for the form y = Acos(x)
  ;X is the independent variable
  ; X = findgen(91)
  ; X = X*!pi/180.
  ;A is an array of free paramenters such that: A[0] = A and A[1] = B
  ;F is the function return, thus F = y.
  ;pder is the derivative of F with respect to each free parameter in A[].
  ; print, 'A = ', A
  ; stop
  F = A[0]*cos(X)
  pder = [cos(X)]
  ; print, 'F = ', F
  stop
end

function fitSine,m0,m90,p0,p90,num,s
  ;stop
  ;m0 is the MEPED 0degree detector electron spectrum
  ;m90 is the same but for the 90degree detector
  ;p0 is the pitch of the 0degree detector
  ;p90 is the pitch of the 90degree detector
  ;stop
  ;flux = [m0,m90]
  ;convert angles to radians
  if p90 gt 90. then tp90 = 90.-(p90-90.) else tp90 = p90 ;pitch needs to be less than 90 and should be reflection in sine curve.
  if p0 gt 90. then tp0 = 90.-(p0-90.) else tp0 = p0
  p90 = tp90
  p0 = tp0

  ;m00 = m0/sin([p0]*(!pi/180.))
  ;m90 = m90/sin([p90]*(!pi/180.))
  if m0 lt 0. then m0 = 0.
  if m90 lt 0. then m90 = 0.

  ; Make sure it fits a sine function. Thus the higher pitch angle needs the higher flux

  if p90 gt p0 and m90 gt m0 then begin

    pitch = [p0, p90]*(!pi/180.)
    ;pitch = reverse(pitch)
    ;pitch = [0.785398, 0.785398]
    flux = [m0, m90]

    A1 = flux[0]/sin(pitch[0]) ;A is array of free parameters, A[0] = A, A[1] = 1
    A2 = flux[1]/sin(pitch[1])
    A = [(A1+A2)/2.]
    ;  stop
    ; A = exp((alog(a1)+alog(a2))/2.)
    ; A = .000001
    ;  stop
    ;print, 'A1 = ', A1
    ;print, 'A2 = ', A2
    ;print, 'Initial A = ', A
    yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fsin',iter=num,status=s)
    A = yfit[0]
    ;yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fcos',iter=num,status=s)
    ;print, 'Final A = ', A
    ;print, ''
    ;print, '90 Degree Detector higher'
    ;stop
  endif

  if p0 gt p90 and m0 gt m90 then begin

    pitch = [p90, p0]*(!pi/180.)
    ;pitch = reverse(pitch)
    ;pitch = [0.785398, 0.785398]
    flux = [m90, m0]
    A1 = flux[1]/sin(pitch[1]) ;A is array of free parameters, A[0] = A, A[1] = 1
    A2 = flux[0]/sin(pitch[0])
    A = [(A1+A2)/2.]
    ;A = .000001
    ;A = exp((alog(a1)+alog(a2))/2.)
    ;print, 'A1 = ', A1
    ;print, 'A2 = ', A2
    ;print, 'Initial A = ', A
    yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fsin',iter=num,status=s)
    A = yfit[0]
    ;yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fcos',iter=num,status=s)
     ; print, 'Final A = ', A
     ; print, yfit
    ;  print, ''
    ; print, '0 Degree detector higher'
;stop
  endif

  if p0 ne p90 and m0 eq m90 then begin

    p0 = p0*(!pi/180.)
    A = m0/sin(p0)
    ;print, 'Fluxes are the same in both detectors'

    if m0 ne 0. then begin

    ;  print, 'Great Scott! The fluxes are the same!
      ;print, m0, m90
      ;stop

    endif

  endif

  if p0 eq p90 and m0 ne m90 then begin

    p0 = p0*(!pi/180.)
    A = m0/sin(p0)
    
    ;print, 'Fluxes are the same in both detectors'

    if m0 ne 0. then begin

   ;   print, 'Great Scott! The pitch angles are the same!
      ;print, m0, m90
      ;stop

    endif

  endif

  if p0 gt p90 and m0 lt m90 then begin

    p0 = p0*(!Pi/180.)
    A = m0/sin(p0)
    ; print, '0 degree angle higher but flux lower'
    ; lower_counter=lower_counter+1

  endif

  if p0 lt p90 and m0 gt m90 then begin

    p0 = p0*(!Pi/180.)
    A = m0/sin(p0)
    ; print, '0 degree angle lower but flux higher'
    ; higher_counter=higher_counter+1

  endif
  
  ;flux = [m0, m90]
  ;pitch = [p0,p90]*(!pi/180.)
  ;w = 1./flux ;using poisson statistics for weighting
  ;w = !null ;using poisson statistics for weighting
  ;fitting to flux = Asin(Bx)
  ;get first guess for A is fit just to m90 and B = 1
  ;A = [flux[1]/sin(pitch[1]),1.0] ;A is array of free parameters, A[0] = A, A[1] = 1
  ;A1 = flux[1]/sin(pitch[1]) ;A is array of free parameters, A[0] = A, A[1] = 1
  ;A2 = flux[0]/sin(pitch[0])
  ;A = [(A1+A2)/2.]
  ;print, 'A1 = ', A1
  ;print, 'A2 = ', A2
  ;print, 'Initial A = ', A
  ;b1 = alog10(A1)
  ;b2 = alog10(A2)
  ;b = [(B1+B2)/2.]
  ;stop
  ;compute best fit
  ; stop
  ;yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fsin',iter=num,status=s)
  ;yfit = CURVEFIT(pitch,flux,w,A,FUNCTION_NAME='fcos',iter=num,status=s)
  ; print, 'Final A = ', A

  ; stop
  ;A_plus_counter = [A, lower_counter, higher_counter]
  return,A
end

pro calculateBLC_v1,indSat

  indSat = 'm02'


  lower_counter=0
  higher_counter=0

  ;this procedure will calculate the bounce loss cone (BLC) of a given satellite (indSat)for each measurement it takes and output it in a new file.
  startSystime = systime() ;for debugging purposes, when did this program start
  sat = 'POES' ;only satellite this program works with, so don't change it.
  ;set up individual satellite data
  if indSat eq 'n15' then begin
    date = ['2003','10','01'] ;year, month, day ;n15
    endDate = ['2003','01','31'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA15'
  endif
  if indSat eq 'n16' then begin
    date = ['2001','01','10'] ;year, month, day ;n16
    endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA16'
  endif
  if indSat eq 'n17' then begin
    date = ['2002','07','12'] ;year, month, day ;n17
    endDate = ['2013','04','10'] ;end date for program ;n17
    fullSatName = 'NOAA17'
  endif
  if indSat eq 'n18' then begin
    date = ['2005','06','07'] ;year, month, day ;n18
    endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA18'
  endif
  if indSat eq 'n19' then begin
    date = ['2009','02','23'] ;year, month, day ;n19
    endDate = ['2009','12','31'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA19'
  endif
  if indSat eq 'm02' then begin
    date = ['2006','12','03'] ;year, month, day ;m02
    endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP02'
  endif
  if indSat eq 'm01' then begin
    date = ['2012','10','03'] ;year, month, day ;m02
    endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP01'
  endif
  if indSat eq 'm03' then begin
    date = ['2019','01','01'] ;year, month, day ;m02
    endDate = ['2019','12','31'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP03'
  endif

  date = ['2014','01','01']   ;starting with where I care.
  enddate = ['2022','12','31']
  fullSatName = 'METOP02'

  fdoydate = JULDAY(fix(date[1]),fix(date[2]),fix(date[0]))
  enddoydate = JULDAY(fix(endDate[1]),fix(endDate[2]),fix(endDate[0]))

  ;loop over days and read in files
  while fdoydate le enddoydate do begin ;loop over days
    print, date
    ;no clobbering
    ;ofilename = 'POES_combinedSpectrum_'+indSat+'_BLC_v2_'+date[0]+date[1]+date[2]+'.nc'
    ;if FILE_TEST('/export/home/pecked/MEPED/data/BLC/'+ofilename) then begin
    ; print, 'This date already processed'
    ; date = advanceDate(fdoydate)
    ; continue
    ;endif
    ;check validity of date

    ;exists1 = FILE_TEST('F:\POES_Data\Level_1_MPE\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc')
    ;exists2 = FILE_TEST('F:\POES_Data\Level_1_MPE\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_90_'+date[0]+date[1]+date[2]+'.nc')
    exists1 = FILE_TEST('D:\POES_Data\MPE_Software\Processed_MEPED_Data\Normalized_nofloor\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc')
    exists2 = FILE_TEST('D:\POES_Data\MPE_Software\Processed_MEPED_Data\Normalized_nofloor\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_90_'+date[0]+date[1]+date[2]+'.nc')


    if exists1 and exists2 then begin

     ; ncid1 = NCDF_OPEN('F:\POES_Data\Level_1_MPE\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc')
     ; ncid2 = NCDF_OPEN('F:\POES_Data\Level_1_MPE\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_90_'+date[0]+date[1]+date[2]+'.nc')
      ncid1 = NCDF_OPEN('D:\POES_Data\MPE_Software\Processed_MEPED_Data\Normalized_nofloor\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc')
      ncid2 = NCDF_OPEN('D:\POES_Data\MPE_Software\Processed_MEPED_Data\Normalized_nofloor\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_90_'+date[0]+date[1]+date[2]+'.nc')

      varid = ncdf_varid(ncid1,'Ecounts')
      ncdf_varget,ncid1,varid,m00
      varid = ncdf_varid(ncid1,'EOcounts')
      ncdf_varget,ncid1,varid,m00_original
      varid = ncdf_varid(ncid1,'Eerror')
      ncdf_varget,ncid1,varid,m00_std
      varid = ncdf_varid(ncid2,'Ecounts')
      ncdf_varget,ncid2,varid,m90
      varid = ncdf_varid(ncid2,'EOcounts')
      ncdf_varget,ncid2,varid,m90_original
      varid = ncdf_varid(ncid2,'Eerror')
      ncdf_varget,ncid2,varid,m90_std
      varid = ncdf_varid(ncid1,'pitch')
      ncdf_varget,ncid1,varid,p00
      varid = ncdf_varid(ncid2,'pitch')
      ncdf_varget,ncid2,varid,p90
      varid = ncdf_varid(ncid1,'foflLat')
      ncdf_varget,ncid1,varid,foflLat
      varid = ncdf_varid(ncid1,'foflLon')
      ncdf_varget,ncid1,varid,foflLon
      varid = ncdf_varid(ncid1,'geogLat')
      ncdf_varget,ncid1,varid,geogLat
      varid = ncdf_varid(ncid1,'geogLon')
      ncdf_varget,ncid1,varid,geogLon
      varid = ncdf_varid(ncid1,'MLT')
      ncdf_varget,ncid1,varid,MLT
      varid = ncdf_varid(ncid1,'lValue')
      ncdf_varget,ncid1,varid,lValue
      varid = ncdf_varid(ncid1,'energy')
      ncdf_varget,ncid1,varid,energy
      varid = ncdf_varid(ncid1,'time')
      ncdf_varget,ncid1,varid,time
      varid = ncdf_varid(ncid1,'rtime')
      ncdf_varget,ncid1,varid,rtime
      varid = ncdf_varid(ncid2,'time')
      ncdf_varget,ncid2,varid,time2
      varid = ncdf_varid(ncid1,'BLC_Angle')
      ncdf_varget,ncid1,varid,BLC_Angle_00
      varid = ncdf_varid(ncid2,'BLC_Angle')
      ncdf_varget,ncid2,varid,BLC_Angle_90
      varid = ncdf_varid(ncid1,'flag')
      ncdf_varget,ncid1,varid,flag1
      varid = ncdf_varid(ncid2,'flag')
      ncdf_varget,ncid2,varid,flag2
      if varid eq -1 then begin
        print, 'No BLC this day'
        NCDF_CLOSE,ncid1
        NCDF_CLOSE,ncid2
        date = advanceDate(fdoydate)
        continue
      endif
      ncdf_varget,ncid1,varid,alpha
      ;stop
      alpha = BLC_Angle_00
      ;a = plot(alpha)

      ;alpha = fltarr(N_ELEMENTS(time))+10. ;random guess to run code with until I have actual BLC_angle. in degrees.
      alpha = alpha*(!pi/180.)
      NCDF_CLOSE,ncid1
      NCDF_CLOSE,ncid2

      ;stop
      ;q = fltarr(2,N_ELEMENTS(time))
      ;blc_flux = dblarr(N_ELEMENTS(energy),N_ELEMENTS(time))
      ;for e = 0,N_ELEMENTS(energy)-1 do begin
      ; for t = 0,N_ELEMENTS(energy)-1 do $
      ;   q[*,t] = fitSine(m00[e,t],m90[e,t],p00[t],p90[t]) ;fit to sine curve
      ; ;calculate BLC with assumptions
      ; blc_flux[e,*] = (q[0,*]*sin(q[1,*]*alpha))*alpha ;not too worried about putting this in the loop since it is small. flux*(alpha-0) = total BLC
      ; stop
      ; ;above assumes that BLC is constant at maximum flux across all angles in the BLC. So just area of a rectangle where one edge is the flux at alpha and the other edge is alpha.
      ;endfor ;energy loop
      ;sm90 = alog(m90) ;log scaled
      ;sm00 = alog(m00)

      if N_ELEMENTS(time) ne N_ELEMENTS(time2) then begin
        print, 'Time oddity has occurred (time discrepancy).'
        ;This happens with extra data that should be taken out.
        if N_ELEMENTS(time) gt N_ELEMENTS(time2) then begin
          print, 'Only using appropriate data from 00 degrees'
          time = time[0:N_ELEMENTS(time2)-1]
          m00 = m00[*,0:N_ELEMENTS(time2)-1]
          m00_original = m00_original[*,0:N_ELEMENTS(time2)-1]
          m00_std = m00_std[*,0:N_ELEMENTS(time2)-1]
          m90 = m90[*,0:N_ELEMENTS(time2)-1]
          m90_original = m90_original[*,0:N_ELEMENTS(time2)-1]
          m90_std = m90_std[*,0:N_ELEMENTS(time2)-1]
          p00 = p00[0:N_ELEMENTS(time2)-1]
          p90 = p90[0:N_ELEMENTS(time2)-1]
          foflLat = foflLat[0:N_ELEMENTS(time2)-1]
          foflLon = foflLon[0:N_ELEMENTS(time2)-1]
          geogLat = geogLat[0:N_ELEMENTS(time2)-1]
          geogLon = geogLon[0:N_ELEMENTS(time2)-1]
          MLT = MLT[0:N_ELEMENTS(time2)-1]
          lValue = lValue[0:N_ELEMENTS(time2)-1]
        endif else begin
          print,'You better get down here boss!' ;this means that 00 has more data than 90
          stop
        endelse
        print, 'Time issue dealt with.'
      endif
      ;stop
      ;-----------------------------------------------
      ; Zero out the spectrum when the original counts are zero
      ;-----------------------------------------------

      for i = 0, N_ELEMENTS(m00[0,*]) -1 do begin

        variable1 = m00_original[0, i]
        variable2 = m90_original[0, i]

        if variable1 eq 0. then begin
          ;stop
          m00[*,i] = 0.

        endif

        if variable2 eq 0. then begin
          ;stop
          m90[*,i] = 0.

        endif

      endfor

      ;stop
      q = fltarr(N_ELEMENTS(time))
      q_total = fltarr(N_ELEMENTS(energy), N_ELEMENTS(time))
      n = fltarr(N_ELEMENTS(time))

      blc_flux = dblarr(N_ELEMENTS(energy),N_ELEMENTS(time))
      blc_std = dblarr(N_ELEMENTS(energy),N_ELEMENTS(time))
      blc_flux_v5 = dblarr(N_ELEMENTS(energy),N_ELEMENTS(time))
      blc_std_v5 = dblarr(N_ELEMENTS(energy),N_ELEMENTS(time))

      for e = 0,N_ELEMENTS(energy)-1 do begin
        for t = 0,N_ELEMENTS(time)-1 do begin
          tm00 = m00[e,t]
          tm90 = m90[e,t]
          tp00 = p00[t]
          tp90 = p90[t]
          ;     stop

          if finite(tm00) eq 0 then tm00=0.
          if finite(tm90) eq 0 then tm90=0.

          if tm00 eq tm90 then begin

            tm90 = tm90+.001
            tp90 = tp90+.001

          endif

          if tm00 gt 1.e10 or tm90 gt 1.e10 then begin
            print, 'Huge Numbers'
            tm00 = 0.
            tm90 = 0.
            tp00 = 0.
            tm90 = 0.
          endif

          ; stop
          if finite(tp90) eq 1 then begin

            q[t] = fitSine(tm00,tm90,tp00,tp90,nt,stat) ;fit to sine curve ;in log space
            q_total[e,t] = q[t]

          endif

          ;  print, q[t], t
          ;stop       ;if stat eq 2 then stop
          ; q[t] = fitSine(m00[e,t],m90[e,t],p00[t],p90[t],nt) ;fit to sine curve ;in log space
          ;     n[t] = nt
          ;     stop
        endfor
        ;        stop
        ; print, "Zero Degree Pitch Angle Higher but Flux lower: " + string(higher_counter)
        ; print, "Zero Degree Pitch Angle Lower but Flux Higher: " + string(lower_counter)

        ;calculate BLC with assumptions
        ;stop


        blc_flux[e,*] = 2.*!pi*q*sin(alpha)/3.
        blc_std[e,*] = SQRT(((m00_std[e,*])^2.) + ((m90_std[e,*])^2.))

        ;blc_flux_v2[e,*] = (q*sin(alpha))*alpha ;not too worried about putting this in the loop since it is small. flux*(alpha-0) = total BLC ;units at this point are c/cm2/s/sr/keV

        ; blc_flux_v3[e,*] = (q*sin(alpha))*!pi ;This is equation 2.3.5 from Evans and Greer.

        ;a = 2.*alpha/!pi ;linear conversion from alpha at satellite to alpha at 120km
        ;blc_flux_v4[e,*] = (q/4.)*((sin((!pi/2.)*(a-2.))/(a-2.))-(sin((!pi/2.)*(a+2.))/(a+2.))) ;integral of Evans and Greer equation 2.3.3 from 0 to 90 degrees.

        ;blc_flux_v5[e,*] = q*SQRT(Bs/B0)*2.*!pi*((sin(!pi/2.))^3.)/3.

        blc_flux_v5[e,*] = 2.*!pi*q*((sin(alpha))^3.)/3. ;integral of equation 2.3.3 from Evans and Greer from 0 to alphaBLC. Correct Method if n=1
        blc_std_v5[e,*] = SQRT(((m00_std[e,*])^2.) + ((m90_std[e,*])^2.))

        ;blc_flux_v6[e,*] = (2.*!pi*q)*(sin(alpha)+((((cos(alpha))^2.)-1.)/2.))  ;above assumes that BLC is constant at maximum flux across all angles in the BLC. So just area of a rectangle where one edge is the flux at alpha and the other edge is alpha.

        ; stop

        if e eq 0 then begin

           ; save, q, alpha, p00, p90, m00, m90, blc_flux, blc_flux_v5, geoglat, geoglon, filename='Math_vars_mod_'+date[0]+date[1]+date[2]+indsat+'.sav'
          ;stop
        endif

        ;stop
      endfor ;energy loop

      ; Build the flag array

      Flag = intarr(N_ELEMENTS(time))

      for ii = 0, N_ELEMENTS(time)-1 do begin

        if flag1[ii] eq 1 or flag2[ii] eq 1 then begin
          Flag[ii] = 1
        endif

      endfor

      ;output nc file
      TimeATTR = 'Time of measurment'
      TimeUnit = 'milliseconds since 1970'
      
      rTimeATTR = 'Time of measurment'
      rTimeUnit = 'hours'

      MLTATTR = 'Satellite Magnetic Local Time'
      MLTUnit = 'hours'

      lValueATTR = 'Satellite L-Value'
      lValueUnit = ' '

      energyATTR = 'Center Energy of Electron/Proton Fluxes'
      energyUnit = 'KeV'

      BLCATTR = 'Bounce Loss Cone Differential Flux'
      BLCUnit = 'counts/cm2/s/keV'

   ;   BLC_stdATTR = 'Bounce Loss Cone Differential Flux 1-sigma Error'
   ;   BLC_stdUnit = 'counts/cm2/s/keV'

      BLCAngleATTR = 'Calculate Bounce Loss Cone Angle'
      BLCAngleUnit = 'degrees'

      geog_latATTR = 'Satellite Geographic Latitude'
      geog_latUnit = 'degrees -90 to 90'

      geog_lonATTR = 'Satellite Geographic Longitude'
      geog_lonUnit = 'degrees 0 to 360'

      fofl_latATTR = 'Satellite Foot of the Field Line Geographic Latitude'
      fofl_latUnit = 'degrees -90 to 90'

      fofl_lonATTR = 'Satellite Foot of the Field Line Geographic Longitude'
      fofl_lonUnit = 'degrees 0 to 360'

      flagATTR = 'Either 0 (no suspicion) or 1 (suspicious) spectrum for 1 or both of the detectors'
      flagUnit = 'unitless'


      ofilename = 'POES_combinedSpectrum_'+indSat+'_BLC_v3_'+date[0]+date[1]+date[2]+'.nc' ;specified ealier.
      ;ofile = NCDF_CREATE('Z:\MEE_Data\Processed_MEPED_Data\Kathy\'+ofilename,/CLOBBER)
      ;ofile = NCDF_CREATE('Z:\MEE_Data\WeiXu_Corrected\'+ofilename,/CLOBBER)
      ;ofile = NCDF_CREATE('F:\POES_Data\Level_2_MPE\'+fullSatName+'\'+ofilename,/CLOBBER)
      ofile = NCDF_CREATE('D:\POES_Data\MPE_Software\Processed_MEPED_Data\Normalized_nofloor\BLC_Flux\'+fullSatName+'\'+ofilename,/CLOBBER)
      ;define dimensions
      energydimID = NCDF_DIMDEF(ofile,'energy',N_ELEMENTS(energy))
      ;edgesdimID = NCDF_DIMDEF(ofile,'edges',N_ELEMENTS(edges))
      timedimID = NCDF_DIMDEF(ofile,'time',/UNLIMITED)
      ;define variables
      timeID = NCDF_VARDEF(ofile,'time',timedimID,/DOUBLE)
      rtimeID = NCDF_VARDEF(ofile,'rtime',timedimID,/DOUBLE)
      geoglatID = NCDF_VARDEF(ofile,'geogLat',timedimID,/DOUBLE)
      geoglonID = NCDF_VARDEF(ofile,'geogLon',timedimID,/DOUBLE)
      fofllatID = NCDF_VARDEF(ofile,'foflLat',timedimID,/DOUBLE)
      fofllonID = NCDF_VARDEF(ofile,'foflLon',timedimID,/DOUBLE)
      MLTID = NCDF_VARDEF(ofile,'MLT',timedimID,/DOUBLE)
      lvalueID = NCDF_VARDEF(ofile,'lValue',timedimID,/DOUBLE)
      energyID = NCDF_VARDEF(ofile,'energy',energydimID,/FLOAT)
      ;edgesID = NCDF_VARDEF(ofile,'edges',edgesdimID,/FLOAT)
      BLCID = NCDF_VARDEF(ofile,'BLC_Flux',[energydimID,timedimID],/DOUBLE)
    ;  BLC_stdID = NCDF_VARDEF(ofile,'BLC_Flux_Error',[energydimID,timedimID],/DOUBLE)
      BLC_AngleID = NCDF_VARDEF(ofile,'BLC_Angle',[timedimID],/double)
      flagID = NCDF_VARDEF(ofile,'Flag',[timedimID],/double)
      ;define attributes
      NCDF_ATTPUT,ofile,timeID,'Full_Name',timeATTR,/CHAR
      NCDF_ATTPUT,ofile,timeID,'Units',timeUnit,/CHAR
      NCDF_ATTPUT,ofile,rtimeID,'Full_Name',rtimeATTR,/CHAR
      NCDF_ATTPUT,ofile,rtimeID,'Units',rtimeUnit,/CHAR
      NCDF_ATTPUT,ofile,geoglatID,'Full_Name',geog_latATTR,/CHAR
      NCDF_ATTPUT,ofile,geoglatID,'Units',geog_latUnit,/CHAR
      NCDF_ATTPUT,ofile,geoglonID,'Full_Name',geog_lonATTR,/CHAR
      NCDF_ATTPUT,ofile,geoglonID,'Units',geog_lonUnit,/CHAR
      NCDF_ATTPUT,ofile,fofllatID,'Full_Name',fofl_latATTR,/CHAR
      NCDF_ATTPUT,ofile,fofllatID,'Units',fofl_latUnit,/CHAR
      NCDF_ATTPUT,ofile,fofllonID,'Full_Name',fofl_lonATTR,/CHAR
      NCDF_ATTPUT,ofile,fofllonID,'Units',fofl_lonUnit,/CHAR
      NCDF_ATTPUT,ofile,MLTID,'Full_Name',MLTATTR,/CHAR
      NCDF_ATTPUT,ofile,MLTID,'Units',MLTUnit,/CHAR
      NCDF_ATTPUT,ofile,lvalueID,'Full_Name',lvalueATTR,/CHAR
      NCDF_ATTPUT,ofile,lvalueID,'Units',lvalueUnit,/CHAR
      NCDF_ATTPUT,ofile,energyID,'Full_Name',energyATTR,/CHAR
      NCDF_ATTPUT,ofile,energyID,'Units',energyUnit,/CHAR
      NCDF_ATTPUT,ofile,BLCID,'Full_Name',BLCATTR,/CHAR
      NCDF_ATTPUT,ofile,BLCID,'Units',BLCUnit,/CHAR
   ;   NCDF_ATTPUT,ofile,BLC_stdID,'Full_Name',BLCATTR,/CHAR
  ;    NCDF_ATTPUT,ofile,BLC_stdID,'Units',BLCUnit,/CHAR
      NCDF_ATTPUT,ofile,BLC_AngleID,'Full_Name',BLCAngleATTR,/CHAR
      NCDF_ATTPUT,ofile,BLC_AngleID,'Units',BLCAngleUnit,/CHAR
      NCDF_ATTPUT,ofile,flagID,'Units',flagUnit,/CHAR
      ;input variables
      NCDF_CONTROL,ofile,/ENDEF
      NCDF_VARPUT,ofile,timeID,time2
      NCDF_VARPUT,ofile,rtimeID,rtime
      NCDF_VARPUT,ofile,geoglatID,geoglat
      NCDF_VARPUT,ofile,geoglonID,geoglon
      NCDF_VARPUT,ofile,fofllatID,fofllat
      NCDF_VARPUT,ofile,fofllonID,fofllon
      NCDF_VARPUT,ofile,MLTID,MLT
      NCDF_VARPUT,ofile,lvalueID,lvalue
      NCDF_VARPUT,ofile,energyID,energy
      NCDF_VARPUT,ofile,BLCID,BLC_flux_v5
   ;   NCDF_VARPUT,ofile,BLC_stdID,BLC_std
      NCDF_VARPUT,ofile,BLC_AngleID,BLC_angle_00
      NCDF_VARPUT,ofile,flagID,flag
      NCDF_CLOSE,ofile
    endif else begin
      if exists1 or exists2 then begin
        print, date
        print, 'Only one angle file for this day'
       ; stop
      endif else print, 'No Data this Day.'
    endelse
    ;advance to the next day
    date = advanceDate(fdoydate)
    ;break ;for debug
  endwhile
  stopSystime = systime() ;for debugging purposes, when did this program stop
  print, 'Start: '+ startSystime
  print, 'Stop: ' + stopSystime

end
