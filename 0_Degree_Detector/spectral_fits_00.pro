;includes
@amoeba_dfpmin_support
@maged_gf
@magpd_gf
@penalty_function_maged
@penalty_function_magpd
@modbess
@POES_spec

;---------------------------------------------------------------------------

; This is the main procedure that begins the process. Set the dates, iroot,
; oroot, satellite, telescope below.  

pro spectral_fits_jmp_v2_00

  common counts, y, dy
  common weighting, dt, energy, Kwfe, Kwfp

  startSystime = systime() ;for debugging purposes, when did this program start
  ;date = ['2019','01','01'] ;year, month, day
  date = ['2019','11','07'] ;year, month, day   flexible start date
  ;date = ['2003','10','20'] ;year, month, day
  ;endDate = ['2003','11','15'] ;end date for program ;n15,16,18,19,m02
  ;date = ['2013','01','01'] ;year, month, day
  ;date = ['1998','07','01'] ;year, month, day ;n15
  ;date = ['2001','01','10'] ;year, month, day ;n16
  ;date = ['2002','07','12'] ;year, month, day ;n17
  ;date = ['2005','06','07'] ;year, month, day ;n18
  ;date = ['2009','02','23'] ;year, month, day ;n19
  ;date = ['2020','01','22'] ;year, month, day ;m02
  ;endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
  endDate = ['2021','12','31'] ;end date for program ;n15,16,18,19,m02
  ;endDate = ['2013','04','10'] ;end date for program ;n17
  ;endDate = ['2014','06','06']    ; flexible end date
  sat = 'POES' ;only satellite this program works with, so don't change it.
  angle = 0 ;look at 0 degree telescopes
  minEnergy = 25. ;minimum energy on the X axis ;I assume this is in keV
  maxElectronEnergy = 10000. ;maximum energy electron retrieval
  maxProtonEnergy = 10000. ;maximum energy proton retrieval
  interpLevel = 27 ;desired resolution of differential particle flux ;number of points in final graph.
  indSat = 'm03'
  fullSatName = 'METOP03'
  
  dy = 0.4
  ;   dy = 0.25
  ;   dt = 60.0
  ;   dt = 1.0
  dt = 16.0 ;using 16 second data for POES ;EDP

  ; MAGED geometrical factors and weighting functions for spectral inversion
  geomfac_ideal = 0.01 ; cm^2 sr ;still true in POES
  ;GEideal = [0.2, 0.5, 1.0, 1.5, 2.5] ; cm^2 sr keV ;will only use bowtie for POES ;EDP
  ;MAGEDEnergyIdeal = [40., 75., 150., 275., 475.] ;same as above comment ;EDP

  ;maged_weighting_function
  ;GEBowtie = [0.212, 0.487, 1.41, 2.72, 4.41]
  ;MAGEDEnergyBowtie = [41.0, 69.1, 141.6, 247.0, 504.7]
  temp = readResponse('Proton','Proton',minEnergy,maxProtonEnergy,maxElectronEnergy,interpLevel)
  ;Kwfe = TRANSPOSE(temp.response)*dt ;units should be cm2*sr*keV*s
  K1 = temp.response*dt ;proton, proton ;cm2*sr*keV*s
  ;Kwfp = temp.response*dt ;units should be cm2*sr*keV*s
  energy = temp.energy
  deltaE = temp.deltaE
  ;stop
  GPBowtie = [0.43,1.35,4.01,11.29,22.03] ;cm2*sr*keV ;going to use protons here (despite this being the variable for MAGED)
  MAGPDEnergyBowtie = [39.,115.,332.,1105.,2723.] ;keV

  temp = readResponse('Proton','Electron',minEnergy,maxProtonEnergy,maxElectronEnergy,interpLevel)
  ;Kwfe = TRANSPOSE(temp.response)*dt ;units should be cm2*sr*keV*s
  K2 = temp.response*dt ;proton, electron
  energy = temp.energy
  deltaE = temp.deltaE

  ;Kwfp = K1+K2 ;units should be cm2*sr*keV*s ;Proton response to electrons and protons
  Kwfp = K1 ;units should be cm2*sr*keV*s ;Proton response to electrons and protons

  temp = readResponse('Electron','Proton',minEnergy,maxProtonEnergy,maxElectronEnergy,interpLevel)
  ;Kwfe = TRANSPOSE(temp.response)*dt ;units should be cm2*sr*keV*s
  K3 = temp.response*dt ;electron, proton
  energy = temp.energy
  deltaE = temp.deltaE

  temp = readResponse('Electron','Electron',minEnergy,maxProtonEnergy,maxElectronEnergy,interpLevel)
  ;Kwfe = TRANSPOSE(temp.response)*dt ;units should be cm2*sr*keV*s
  K4 = temp.response*dt ;electron, electron
  energy = temp.energy
  deltaE = temp.deltaE


  ; MAGPD geometrical factors and weighting functions for spectral inversion
  ;GPideal = [0.3, 0.6, 0.8, 1.0, 4.5] ;only using bowtie ;EDP
  ;MAGPDEnergyIdeal = [95., 140., 210., 300., 575.] ;only using bowtie ;EDP

  ;MAGPDmeasured = 1

  ;if MAGPDmeasured eq 0 then begin
  ; bowtie results from recommended GFs
  ;   magpd_weighting_function, MEASURED = 0
  ;   GPBowtie = [0.242, 0.526, 0.724, 0.925, 3.96]
  ;   MAGPDEnergyBowtie = [92.7, 133.2, 201.5, 290.4, 471.8]
  ;endif else begin
  ; bowtie results from measured GFs
  ;   magpd_weighting_function, MEASURED = 1
  ;   GPBowtie = [0.163, 0.519, 0.843, 0.838, 3.71]
  ;   MAGPDEnergyBowtie = [94.4, 137.8, 200.2, 288.4, 470.0]
  ;endelse
  ;POES_electron_weighting_function
  GEBowtie = [1.12,2.26,2.44,3.58] ;these are the MEPED electron channels E1-E4.
  ;GEBowtie = [1.12,2.26,2.44] ;these are the MEPED electron channels E1-E3.
  MAGEDEnergyBowtie = [72.,193.,419.,879.]
  ;MAGEDEnergyBowtie = [72.,193.,419.]


  ;Create loop of days
  fdoydate = JULDAY(fix(date[1]),fix(date[2]),fix(date[0]))
  enddoydate = JULDAY(fix(endDate[1]),fix(endDate[2]),fix(endDate[0]))
  endOfYear = 1
  noClobber = 1
  oneDay = 0 ;if oneDay then one day is run only instead of a full year.

  while fdoydate le enddoydate do begin

    Wpl_total = []
    wRM_total = []
    wDM_total = []
    wEE_total = []

   print, date
    
 ;--------------------------------------------------------------------------------------------
 ; Read in count data from MEPED   --- As of 2015 the data format changed
 ;--------------------------------------------------------------------------------------------
    
    if date[0] ge 2012 then begin
      
      temp = pCounts2(angle,'Electron',sat,date) ; reads new netcdf files
      if (temp.particleCount)[0] eq -999 then begin
          print, 'Bad data this day.'
          date = advanceDate(fdoydate)
          continue
          stop
      endif
      
    ENDIF ELSE BEGIN
        
      temp = pCounts(angle,'Electron',sat,date) ;read CDF file
      if (temp.particleCount)[0] eq -999 then begin
          print, 'Bad data this day.'
          date = advanceDate(fdoydate)
          continue
      endif
      
    endelse
    
;    stop
    electronCount = temp.particleCount
    geogLL = temp.geogLL ;geographic latitude and longitude of the satellite. Not really helpful.
    time = temp.time
    lValue = temp.lValue ;can be converted into invariant latitude for estimate of magnetic latitude.
    foflLL = temp.foflLL ;2 dimensions, first dimension declares lat [0,*] or lon [1,*]. Also true for geogLL. Foot of the field line geographic lat and lon. This is where the particles are coming into the atmosphere in geographic coordinates.
    MLT = temp.MLT ;magnetic local time of the satellite.
    pitch = temp.pitch ;pitch angle of measured particles
    bfofl = temp.bfofl ; magnetic field strength at foot of satellite
    blocal = temp.blocal ; magnetic field strength at satellite
    blc_angle = temp.blc_angle ; calculated bounce loss cone angle
    ;protonCount = temp.particleCount
    restore, 'D:\POES_Data\MPE_Software\0_Degree_Detector\Pcounts_temp.sav'
    ;electronCount = particlecount
    protonCount = pparticlecount
    geogLL = temp.geogLL
    time = temp.time
    msec_time = msec_time
    ;    stop
    
    ;oldyr = fix(date[0])
    ;date = advanceDate(fdoydate)
    ;fdoydate = fdoydate+1.
    ;CALDAT,fdoydate,dmo,ddy,dyr
    ;date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]

    ;if oldyr ne dyr then endOfYear = 0 else date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]
    ;stop

    ;we will first use protons, P1-P5. We after solution is decided on, we will use that to calculate a new P6. Once that is done, we can subtract it from the original P6 to find out what the E4 channel is.
    Nstart = 0 ;will always start with the first time point.
    Nend = N_ELEMENTS(time) ;end with the last time
    N1min = N_ELEMENTS(time) ;time variable
    ; stop
    ;   MAGEDMaxwFirst = fltarr(9, N1min, 2)
    ;dimensions are [Number of telescopes, time, Nq]. So if using only P1-P5 (Ny = 5 channels and Nq = 2 parameterizations (q)), then dimensions are [1,time,2]. ;EDP
    MAGPDRMQ = fltarr(N1min, 2)
    MAGPDRMCovQQ = fltarr(N1min, 2, 2)
    MAGPDDMQ = fltarr(N1min, 4)
    MAGPDDMCovQQ = fltarr(N1min, 4, 4)
    MAGPDPLQ = fltarr(N1min, 2)
    MAGPDPLCovQQ = fltarr(N1min, 2, 2)
    MAGPDEEQ = fltarr(N1min, 2)
    MAGPDEECovQQ = fltarr(N1min, 2, 2)
    jpflux = fltarr(5, N1min) ;best guess fluxes using G from bowtie ;EDP ;this has dimensinos [Ny,time] ;I probably don't need to track jflux since it is the first guess.
    wkp = fltarr(N1min, 4) ;this has dimensions [time, number of telescopes, number of spectral functions].

    for i = Nstart, Nend-1 do begin ;time loop ;EDP
     ;for i = 4410, 4411 do begin ;time loop ;EDP
      ;print, 'Nstart: '+string(Nstart)+' Nend: '+string(Nend)+' i: '+string(i)
      ;print, "Proton Number :", i
      ;--------------------------------------------------------------------------------------------
      ; Proton differential flux calculation
      ;--------------------------------------------------------------------------------------------

      jpflux[*,i] = protoncount[0:4,i]/(GPbowtie) ;counts/sec/cm2/sr/keV ; y = Kx, x = jflux, K = G = GPbowtie (cm2srkev)
      y = protoncount[0:4,i]*dt ;counts/sec * sec = counts
      ; relativistic Maxwellian
      Q0RM = q_first_guess_RMp(jpflux[*,i], MAGPDEnergyBowtie)
      QsRM = amoeba(0.001, FUNCTION_NAME = 'penalty_magpd_relmaxwell', $
        P0 = Q0RM, SCALE = [0.2*abs(Q0RM[0]), 0.5*abs(Q0RM[1])])
      if n_elements(QsRM) eq 1 then begin
        ;print, i, QsRM
        ;PRINT, y
        covqqRM = [[-1.0e31, -1.0e31], [-1.0e31, -1.0e31]]
        ;stop
      endif else begin
        ell_RM = penalty_function_magpd(QsRM, 'relmaxwellp')
        dldq_RM = penalty_function_magpd_gradient(QsRM, $
          'relmaxwellp', 'relmaxwellp_gradient' )
        d2ldq2RM = penalty_function_magpd_hessian(QsRM, $
          'relmaxwellp', 'relmaxwellp_gradient', 'relmaxwellp_curvature')
        covqqRM = covariance_qq(d2ldq2RM)
      endelse
      ;   stop
      ; power law
      Q0PL = q_first_guess_PL(jpflux[*,i], MAGPDEnergyBowtie)
      QsPL = amoeba(0.001, FUNCTION_NAME = 'penalty_magpd_powerlaw', $
        P0 = Q0PL, SCALE = [0.2*abs(Q0PL[0]), 0.5*abs(Q0PL[1])])
      if n_elements(QsPL) eq 1 then begin
        ;print, i, QsPL
        ;PRINT, y
        covqqPL = [[-1.0e31, -1.0e31], [-1.0e31, -1.0e31]]
      endif else begin
        ell_PL = penalty_function_magpd(QsPL, 'powerlaw')
        dldq_PL = penalty_function_magpd_gradient(QsPL, $
          'powerlaw', 'powerlaw_gradient' )
        d2ldq2PL = penalty_function_magpd_hessian(QsPL, $
          'powerlaw', 'powerlaw_gradient', 'powerlaw_curvature')
        covqqPL = covariance_qq(d2ldq2PL)
      endelse
      ;   stop
      ; energy exponential
      Q0EE = q_first_guess_EE(jpflux[*,i], MAGPDEnergyBowtie)
      QsEE = amoeba(0.001, FUNCTION_NAME = 'penalty_magpd_exponential', $
        P0 = Q0EE, SCALE = [0.2*abs(Q0EE[0]), 0.5*abs(Q0EE[1])])
      if n_elements(QsEE) eq 1 then begin
        ;print, i, QsEE
        ;PRINT, y
        covqqEE = [[-1.0e31, -1.0e31], [-1.0e31, -1.0e31]]
      endif else begin
        ell_EE = penalty_function_magpd(QsEE, 'energyexponential')
        dldq_EE = penalty_function_magpd_gradient(QsEE, $
          'energyexponential', 'energyexponential_gradient' )
        d2ldq2EE = penalty_function_magpd_hessian(QsEE, $
          'energyexponential', 'energyexponential_gradient', 'energyexponential_curvature')
        covqqEE = covariance_qq(d2ldq2EE)
      endelse
      ;stop
      ; Double Relativistic Maxwellian Spectrum
      Q0DM = q_first_guess_DMp(jpflux[*,i], MAGPDEnergyBowtie)
      QsDM = amoeba(0.001, FUNCTION_NAME = 'penalty_magpd_drelmaxwell', $
        P0 = Q0DM, SCALE = [0.2*abs(Q0DM[0]), 0.5*abs(Q0DM[1]), 0.2*abs(Q0DM[2]), 0.5*abs(Q0DM[3])])
      if n_elements(QsDM) eq 1 then begin
        ;print, i, QsDM
        ;PRINT, y
        covqqDM = [[-1.0e31, -1.0e31], [-1.0e31, -1.0e31]]
      endif else begin
        ell_DM = penalty_function_magpd(QsDM, 'fdrelmaxwellp')
        dldq_DM = penalty_function_magpd_gradient(QsDM, $
          'fdrelmaxwellp', 'fdrelmaxwellp_gradient' )
        d2ldq2DM = penalty_function_magpd_hessian(QsDM, $
          'fdrelmaxwellp', 'fdrelmaxwellp_gradient', 'fdrelmaxwellp_curvature')
        covqqDM = covariance_qq(d2ldq2DM)
      endelse
      ;   stop
      if n_elements(QsRM) gt 1 then begin
        MAGPDRMQ[i, *] = QsRM
        MAGPDRMCovQQ[i, *, *] = covqqRM
        wRM = exp(-ell_RM[5,0]-2.0)
      endif else begin
        MAGPDRMQ[i, *] = [-1.0e31, -1.0e31]
        MAGPDRMCovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wRM = 0.0
      endelse
      if n_elements(QsDM) gt 1 then begin
        MAGPDDMQ[i, *] = QsDM
        MAGPDDMCovQQ[i, *, *] = covqqDM
        ;wDM = exp(-ell_DM[5,0]-2.0)
        wDM = exp(-ell_DM[5,0]-4.0)
      endif else begin
        MAGPDDMQ[i, *] = dblarr(4)-1.0e31
        MAGPDDMCovQQ[i, *, *] = dblarr(4,4)-1.0e31
        wDM = 0.0
      endelse
      ;         MAGPDNcalls[j, i] = MENcalls
      if n_elements(QsPL) gt 1 then begin
        MAGPDPLQ[i, *] = QsPL
        MAGPDPLCovQQ[i, *, *] = covqqPL
        wPL = exp(-ell_PL[5,0]-2.0)
      endif else begin
        MAGPDPLQ[i, *] = [-1.0e31, -1.0e31]
        MAGPDPLCovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wPL = 0.0
      endelse

      if n_elements(QsEE) gt 1 then begin
        MAGPDEEQ[i, *] = QsEE
        MAGPDEECovQQ[i, *, *] = covqqEE
        wEE = exp(-ell_EE[5,0]-2.0)
      endif else begin
        MAGPDEEQ[i, *] = [-1.0e31, -1.0e31]
        MAGPDEECovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wEE = 0.0
      endelse
      wkp[i,*] = [wRM, wPL, wEE, wDM]/(wRM + wPL + wEE+wDM)
      ;stop
    endfor


    print, "Beginning plot_fit_pflux"
    plot_fit_pflux, Nend, MAGPDEnergyBowtie, GPBowtie, jpflux, MAGPDEnergyBowtie, MAGPDRMQ, MAGPDRMCovQQ, $
      MAGPDPLQ, MAGPDPLCovQQ, MAGPDEEQ, MAGPDEECovQQ, MAGPDDMQ, MAGPDDMCovQQ, wkp, protoncount,protonFluxes, pLambda, pSigma
    ;stop
    print, "plot_fit_pflux complete"
    ;--------------------------------------------------------------------------------------------
    ; Proton Contamination Removal
    ;--------------------------------------------------------------------------------------------
    
    ;Get proton contamination counts in electron channels and subtract from E1-E3.
    protoncontamination = dblarr(3,Nend)
    new_electronCount = dblarr(3,Nend) ;E1,E2,E3 x time
    for i = Nstart, Nend-1 do begin ;time loop ;EDP
      pContam = K3##protonFluxes[*,i] ;should return proton contamination
      protoncontamination[*,i] = reform(pcontam/dt)
      new_electronCount[*,i] = electronCount[*,i] - REFORM(pContam/dt)
    endfor
    
    ;--------------------------------------------------------------------------------------------
    ; E4 Channel added
    ;--------------------------------------------------------------------------------------------
    
    ;electron counts is holding all the electron data for E1-E3. Need to add in E4.
    new_electronCount = [new_electronCount,(protonCount[5,*]-(pLambda[5,*]/dt))] ;P6 - calculated P6 = E4
    ;stop ;check for negative cases. I may need to add in a negativity clause of where lt 0 set to 0.

    badData = where(new_electroncount lt 0.)
    new_electronCount[badData] = 0.
    Kwfe = [[K4],[K2[*,5]]] ;proton telescope response to electrons in P6 is the response function for E4.
 
 
count_rate_issues = 'yes'

if count_rate_issues eq 'yes' then begin    
    ;------------------------------------------------------------------------------------
    ; Deal with count rate channel issues
    ;------------------------------------------------------------------------------------

    for index = 0, N_ELEMENTS(new_electroncount[0, *])-1 do begin
      if new_electroncount[0, index] gt 2.5 then begin
        if new_electroncount[2, index] eq 0. or new_electroncount[2, index] lt new_electroncount[3, index] then begin
          if new_electroncount[1, index] ne 0. then begin;
          ;stop
            ;--------
            ; If E3 is 0 and E4 is not zero, take logarthmic mean of E2 and E4 for E3
            ;--------
            if new_electroncount[3, index] ne 0. then begin
 ;             ;print, electroncount[*, index]
 ;             ;print, new_electroncount[*, index]
 ;             ;new_electroncount[2, index] = exp((alog(new_electroncount[1, index])+alog(new_electroncount[3, index]))/2.)
              new_electroncount[*, index] = 0.
 ;             ;print, new_electroncount[*, index]
          ;stop
            endif
          endif
        endif
      endif

    endfor
    
    ;------------------------------------------------------------------------------------
    ; Deal with count rate channel issues
    ;------------------------------------------------------------------------------------
  ;  counter = 0
    for index = 0, N_ELEMENTS(new_electroncount[0, *])-1 do begin
      if new_electroncount[0, index] gt 2.5 then begin
        if new_electroncount[1, index] lt new_electroncount[3, index] then begin

          ;--------
          ; If E2 is lt and E4 zero channels
          ;--------

          ;print, electroncount[*, index]
          ;print, new_electroncount[*, index]
          new_electroncount[*, index] = 0.
          ;print, new_electroncount[*, index]
          ;stop
   ;       counter = counter + 1

        endif
      endif
    endfor
    
    ;------------------------------------------------------------------------------------
    ; Deal with count rate channel issues
    ;------------------------------------------------------------------------------------
    ;  counter = 0
    for index = 0, N_ELEMENTS(new_electroncount[0, *])-1 do begin
      if new_electroncount[0, index] gt 2.5 then begin
        if new_electroncount[0, index] lt new_electroncount[3, index] then begin

          ;--------
          ; If E1 is lt and E4 zero channels
          ;--------

          ;print, electroncount[*, index]
          ;print, new_electroncount[*, index]
          new_electroncount[*, index] = 0.
          ;print, new_electroncount[*, index]
          ;stop
          ;       counter = counter + 1

        endif
      endif
    endfor      

endif        

   ; stop
    ;------------------------------------------------------------------------------------
    ; Build a noise floor by removing any of the data with less than .125 counts/sec ;JMP
    ;------------------------------------------------------------------------------------
print, 'Noise_Floor turned off'
    ;------------------------------------------------------------------------------------
    ;print, N_ELEMENTS(where(new_electroncount[0, *] eq 0.))
    ;NOISE_FLOOR_INDEX = where(new_electroncount[0, *] le 2.5)
    ;new_electroncount[*, NOISE_FLOOR_INDEX] = 0.
    ;print, N_ELEMENTS(where(new_electroncount[0, *] eq 0.))
    ;------------------------------------------------------------------------------------


 ; stop
    ;new_electronCount = electronCount
    ;Kwfe = K4
    ;useE3 = 1
    ;useE4 = 0
    ;if useE3 eq 1 and useE4 eq 0 then begin
    ;  Kwfe = K4 ;if we are ignoring E4, then we do not use the E4 weighting function.
    ;  MAGEDEnergyBowtie = MAGEDEnergyBowtie[0:2]
    ;endif
    ;if useE4 eq 1 and useE3 eq 0 then begin
    ;  Kwfe = [[K4[*,0:1]],[K2[*,5]]] ;if we are ignoring E3, then we do not use the E3 weighting function
    ;  MAGEDEnergyBowtie = [MAGEDEnergyBowtie[0:1],MAGEDEnergyBowtie[3]]
    ;endif
    ;if useE3 eq 0 and useE4 eq 0 then begin
    ;  print, 'Either E3 or E4 must be turned on, else this program will crash. Choose one.'
    ;  stop
    ;endif
    ;if useE3 eq 1 and useE4 eq 1 then begin
    ;  print, 'E3 and E4 cannot both be used as they bicker too much. Choose one.'
    ;  stop
    ;endif

    ;stop

    ;Break E1-E4 into spectra.
    MAGEDRMQ = fltarr(N1min, 2)
    MAGEDRMCovQQ = fltarr(N1min, 2, 2)
    MAGEDPLQ = fltarr(N1min, 2)
    MAGEDPLCovQQ = fltarr(N1min, 2, 2)
    MAGEDEEQ = fltarr(N1min, 2)
    MAGEDEECovQQ = fltarr(N1min, 2, 2)
    MAGEDDMQ = fltarr(N1min, 4)
    MAGEDDMCovQQ = fltarr(N1min, 4, 4)
    jflux = fltarr(4, N1min)
    wk = fltarr(N1min, 4)
    ;jflux = fltarr(3, N1min) ;assuming only using 3 channels, not four.
    ;wk = fltarr(N1min, 3)

    for i = Nstart, Nend-1 do begin
      
      
      ;print, "Electron Number :", i
      ;for i = 4410, 4411 do begin
      ;for i = 1, 100 do begin
      ;print, 'Nstart: '+string(Nstart)+' Nend: '+string(Nend)+' i: '+string(i)
      ;if useE3 eq 1 and useE4 eq 0 then begin
      ;  jflux[*,i] = new_electronCount[0:2,i]/(GEbowtie[0:2])
      ;  y = new_electronCount[0:2,i]*dt
      ;endif
      ;if useE4 eq 1 and useE3 eq 0 then begin
      ;  jflux[*,i] = [new_electronCount[0:1,i],new_electronCount[3,i]]/([GEbowtie[0:1],GEbowtie[3]])
      ;  y = [new_electronCount[0:1,i],new_electronCount[3,i]]*dt
      ;endif
      ;if useE3 eq 1 and useE4 eq 1 then begin
      ;  print, 'I already made it clear that E3 and E4 do not get along. Now the program will crash. Thanks a lot!'
      ;  stop

      ;--------------------------------------------------------------------------------------------
      ; Electron differential flux calculation
      ;--------------------------------------------------------------------------------------------

      jflux[*,i] = new_electronCount[*,i]/(GEbowtie) ;counts/sec/cm2/sr/keV ; y = Kx, x = jflux, K = G = GPbowtie (cm2srkev)
      y = new_electronCount[*,i]*dt ;counts/sec * sec = counts
      chan_energy = [72.000, 193.000, 419.000, 879.000]

      ;endif
      ;stop
      ; relativistic Maxwellian
      Q0RM = q_first_guess_RM(jflux[*,i], MAGEDEnergyBowtie)
      QsRM = amoeba(0.001, FUNCTION_NAME = 'penalty_maged_relmaxwell', $
        P0 = Q0RM, SCALE = [0.2*abs(Q0RM[0]), 0.5*abs(Q0RM[1])])
      if n_elements(QsRM) gt 1 then begin
        ell_RM = penalty_function_maged(QsRM, 'relmaxwell');,1)
        dldq_RM = penalty_function_maged_gradient(QsRM, $
          'relmaxwell', 'relmaxwell_gradient' )
        d2ldq2RM = penalty_function_maged_hessian(QsRM, $
          'relmaxwell', 'relmaxwell_gradient', 'relmaxwell_curvature')
        covqqRM = covariance_qq(d2ldq2RM)
        MAGEDRMQ[i, *] = QsRM
        MAGEDRMCovQQ[i, *, *] = covqqRM
        wRM = exp(-ell_RM[4,0]-2.0)
      endif else begin
        MAGEDRMQ[i-Nstart, *] = [-1.0e31, -1.0e31]
        MAGEDRMCovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wRM = 0.0
      endelse
      ;stop
      ; power law
      Q0PL = q_first_guess_PL(jflux[*,i], MAGEDEnergyBowtie)
 ;     stop
      QsPL = amoeba(0.001, FUNCTION_NAME = 'penalty_maged_powerlaw', $
        P0 = Q0PL, SCALE = [0.2*abs(Q0PL[0]), 0.5*abs(Q0PL[1])]);,NCALLS=ncalls, NMAX = long(30000))
 ;     stop
      if n_elements(QsPL) gt 1 then begin
        ell_PL = penalty_function_maged(QsPL, 'powerlaw')
        dldq_PL = penalty_function_maged_gradient(QsPL, $
          'powerlaw', 'powerlaw_gradient' )
        d2ldq2PL = penalty_function_maged_hessian(QsPL, $
          'powerlaw', 'powerlaw_gradient', 'powerlaw_curvature')
        covqqPL = covariance_qq(d2ldq2PL)
        MAGEDPLQ[i, *] = QsPL
        MAGEDPLCovQQ[i, *, *] = covqqPL
        wPL = exp(-ell_PL[4,0]-2.0)
      endif else begin
        MAGEDPLQ[i, *] = [-1.0e31, -1.0e31]
        MAGEDPLCovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wPL = 0.0
      endelse
      ;stop
      ; energy exponential
      Q0EE = q_first_guess_EE(jflux[*,i], MAGEDEnergyBowtie)
      QsEE = amoeba(0.001, FUNCTION_NAME = 'penalty_maged_exponential', $
        P0 = Q0EE, SCALE = [0.2*abs(Q0EE[0]), 0.5*abs(Q0EE[1])])
      if n_elements(QsEE) gt 1 then begin
        ell_EE = penalty_function_maged(QsEE, 'energyexponential')
        dldq_EE = penalty_function_maged_gradient(QsEE, $
          'energyexponential', 'energyexponential_gradient' )
        d2ldq2EE = penalty_function_maged_hessian(QsEE, $
          'energyexponential', 'energyexponential_gradient', 'energyexponential_curvature')
        covqqEE = covariance_qq(d2ldq2EE)
        MAGEDEEQ[i, *] = QsEE
        MAGEDEECovQQ[i, *, *] = covqqEE
        wEE = exp(-ell_EE[4,0]-2.0)
      endif else begin
        MAGEDEEQ[i, *] = [-1.0e31, -1.0e31]
        MAGEDEECovQQ[i, *, *] = [[-1.0e31, -1.0e31],[-1.0e31, -1.0e31]]
        wEE = 0.0
      endelse
      ;stop
      ; Double Relativistic Maxwellian Spectrum
      Q0DM = q_first_guess_DM(jflux[*,i], MAGPDEnergyBowtie)
      QsDM = amoeba(0.001, FUNCTION_NAME = 'penalty_maged_drelmaxwell', $
        P0 = Q0DM, SCALE = [0.2*abs(Q0DM[0]), 0.5*abs(Q0DM[1]), 0.2*abs(Q0DM[2]), 0.5*abs(Q0DM[3])])
      if n_elements(QsDM) eq 1 then begin
        ;print, i, QsDM
        ;PRINT, y
        covqqDM = [[-1.0e31, -1.0e31], [-1.0e31, -1.0e31]]
        MAGEDDMQ[i,*] = [-1.0e31, -1.0e31, -1.0e31, -1.0e31]
        MAGEDDMCovQQ[i,*,*] = [[-1.0e31, -1.0e31, -1.0e31, -1.0e31], [-1.0e31, -1.0e31, -1.0e31, -1.0e31], [-1.0e31, -1.0e31, -1.0e31, -1.0e31], [-1.0e31, -1.0e31, -1.0e31, -1.0e31]]
      endif else begin
        ell_DM = penalty_function_maged(QsDM, 'fdrelmaxwell')
        dldq_DM = penalty_function_maged_gradient(QsDM, $
          'fdrelmaxwell', 'fdrelmaxwell_gradient' )
        d2ldq2DM = penalty_function_maged_hessian(QsDM, $
          'fdrelmaxwell', 'fdrelmaxwell_gradient', 'fdrelmaxwell_curvature')
        covqqDM = covariance_qq(d2ldq2DM)
        MAGEDDMQ[i,*] = QsDM
        MAGEDDMCovQQ[i,*,*] = covqqDM
        ;wDM = exp(-ell_DM[3,0]-2.0)
        wDM = exp(-ell_DM[4,0]-4.0) ;4.0 for the 4 free parameters, versus 2.0 for 2 free parameters.
      endelse

      wRM_total = [wRM_total, wRM]
      wDM_total = [wDM_total, wDM]
      wPL_total = [WPL_total, wPL]
      wEE_total = [wEE_total, wEE]

      ; calculate weighted combination
      wk[i,*] = [wRM, wPL, wEE, wDM]/(wRM + wPL + wEE + wDM)
     
      ;stop
      
    endfor
;stop
  print, "Beginning plot_fit_eflux"
    plot_fit_eflux, Nend, MAGEDEnergyBowtie, GEBowtie, jflux, MAGEDRMQ, MAGEDRMCovQQ, $
      MAGEDPLQ, MAGEDPLCovQQ, MAGEDEEQ, MAGEDEECovQQ, MAGEDDMQ, MAGEDDMCovQQ, wk, new_electronCount, $
      electronFluxes, eLambda, eLambda_pl,elambda_RM,elambda_DM,elambda_EE, eSigma
  print, "plot_Fit_eflux Complete"

;stop
    ;stop

    ;new_electronCount[badData] = -999.
    protonCount[where(protonCount lt 0)] = -999.
    electronCount[where(electronCount lt 0)] = -999.
    protonFluxes[where(protonFluxes lt 0)] = -999.
    electronFluxes[where(electronFluxes lt 0)] = -999.
    eLambda[where(eLambda lt 0)] = -999.
    pLambda[where(pLambda lt 0)] = -999.
    eSigma[where(eSigma lt 0)] = -999.
    pSigma[where(pSigma lt 0)] = -999.
    eSigma[where(finite(eSigma) eq 0)] = -999. ;bad data gives a NaN for sigma, so turning those spots into -999.
    pSigma[where(finite(pSigma) eq 0)] = -999.

    ; stop
    
    ;--------------------------------------------------------------------------------------------
    ; Output Flux Data
    ;--------------------------------------------------------------------------------------------
     print, "Saving File"
    ;outputNCfile,energy,electronFluxes,protonFluxes,eSigma,pSigma,eLambda,pLambda,geogLL,foflLL,time,MLT,lValue,electronCount,protonCount,new_electronCount,MAGEDRMQ,MAGEDPLQ,MAGEDEEQ,MAGEDDMQ,MAGPDRMQ,MAGPDPLQ,MAGPDEEQ,MAGPDDMQ,'m02','METOP02',date,dt,angle
    ;outputNCfile,energy,electronFluxes,protonFluxes,eSigma,pSigma,eLambda,pLambda,geogLL,foflLL,time,MLT,lValue,electronCount,protonCount,new_electronCount,MAGEDRMQ,MAGEDPLQ,MAGEDEEQ,MAGEDDMQ,MAGPDRMQ,MAGPDPLQ,MAGPDEEQ,MAGPDDMQ,'n17','METOP02',date,dt,angle
    ;outputNCfile,energy,electronFluxes,protonFluxes,eSigma,pSigma,eLambda,pLambda,geogLL,foflLL,time,MLT,lValue,pitch,electronCount,protonCount,new_electronCount,MAGEDRMQ,MAGEDPLQ,MAGEDEEQ,MAGEDDMQ,MAGPDRMQ,MAGPDPLQ,MAGPDEEQ,MAGPDDMQ,wk,wkp,indSat,fullSatName,date,dt,angle,msec_time
    ;outputNCfile,energy,electronFluxes,protonFluxes,eSigma,pSigma,eLambda,pLambda,geogLL,foflLL,time,MLT,lValue,pitch,electronCount,protonCount,new_electronCount,MAGEDRMQ,MAGEDPLQ,MAGEDEEQ,MAGEDDMQ,MAGPDRMQ,MAGPDPLQ,MAGPDEEQ,MAGPDDMQ,wk,wkp,indSat,fullSatName,date,dt,angle              
    outputNCfile,energy,electronFluxes,protonFluxes,eSigma,pSigma,eLambda,eLambda_pl,eLambda_dm,eLambda_rm,eLambda_ee,pLambda,geogLL,foflLL,time,MLT,lValue,pitch,electronCount,protonCount,new_electronCount,MAGEDRMQ,MAGEDPLQ,MAGEDEEQ,MAGEDDMQ,MAGPDRMQ,MAGPDPLQ,MAGPDEEQ,MAGPDDMQ,wk,wkp,indSat,fullSatName,date,dt,angle,msec_time,protoncontamination,bfofl,blocal,blc_angle

    if oneDay then break ;should break out of while loop
    ;stop
    ;advance day
    ;oldyr = fix(date[0])
    date = advanceDate(fdoydate)
    ;fdoydate = fdoydate+1.
    ;CALDAT,fdoydate,dmo,ddy,dyr
    ;;if oldyr ne dyr then endOfYear = 0 else date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]
    ;date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]
   ; stop
  endwhile
  stopSystime = systime() ;for debugging purposes, when did this program stop
  print, 'Start: '+ startSystime
  print, 'Stop: ' + stopSystime
end



PRO ColorCat12Png, ColorTrue

; Categorical color key from U. Oregon Datagraphics
; http://geography.uoregon.edu/datagraphics/color/Cat_12.txt
; See Light and Bartlein, Eos, v. 85, no. 40, p. 385, October 5, 2004.
; Produces TrueColor colors in 12-element vector

RGB256 = LONARR(3, 13)

RGB256[*,  0] = [  0,   0,   0] ; black
RGB256[*,  1] = [255, 191, 127] ; light orange
RGB256[*,  2] = [255, 127,   0] ; dark orange
RGB256[*,  3] = [255, 255, 153] ; light yellow
RGB256[*,  4] = [255, 255,  50] ; dark yellow
RGB256[*,  5] = [178, 255, 140] ; light green
RGB256[*,  6] = [ 50, 255,   0] ; dark green
RGB256[*,  7] = [165, 237, 255] ; light blue
RGB256[*,  8] = [ 25, 178, 255] ; dark blue
RGB256[*,  9] = [204, 191, 255] ; light purple
RGB256[*, 10] = [101,  76, 255] ; dark purple
RGB256[*, 11] = [255, 153, 191] ; light red
RGB256[*, 12] = [229,  25,  50] ; dark red

ColorTrue = RGB256[0,*] + 256L * (RGB256[1,*] + 256L * RGB256[2,*])

!P.BACKGROUND = 2^24-1
!P.COLOR = 0

END


;----------------------------------------------------------------------------------

function create_jday_1min, yr, mo, dy, start_hr, end_hr

; create julian dates at 1 min intervals centered on 30 sec

   Nmin = (end_hr-start_hr)*60L
   jday_1min = dblarr(Nmin)
   for hr = start_hr, end_hr-1 do begin
      for min = 0, 59 do begin
         jday_1min[(hr-start_hr)*60L + min] = julday(mo, dy, yr, hr, min, 30)
      endfor
   endfor
   
   return, jday_1min

end


;----------------------------------------------------------------------------------

FUNCTION xticksMAGED, axis, index, value

caldat, value, month,day,year,hour,minute,sec

return, string(hour,minute,format="(i2.2,':',i2.2)")

end


;----------------------------------------------------------------------------------

pro plot_flux_mag, Flux1, magJday, magArray, plot_type, StartJul, EndJul, Yrng, $
   nxticks, nxminor, out_fn_core, DateStr, flux1_title, mag_title, $
   FLUX2STRUCT = Flux2, YRNG2 = Yrng2, FLUX2TITLE = flux2_title

   if plot_type eq 'png' then begin
      set_plot, 'win'
      if keyword_set(Flux2) then begin
         WINDOW, 20, XPOS = 0, YPOS = 0, XSize = 600, YSize = 700
         !P.MULTI=[0,1,3] 
         ChSz = 2.0
      endif else begin
         WINDOW, 20, XPOS = 0, YPOS = 0, XSize = 600, YSize = 500
         !P.MULTI=[0,1,2]  
         ChSz = 1.0     
      endelse
      !P.charthick = 1
      !P.thick = 1    
      ColorCat12Png, color_scale
      SatColor = [color_scale[0], color_scale[12], color_scale[8]]
      !X.ticklen = 0.08
      !Y.ticklen = 0.02   
   endif
   
   caldat, StartJul, Sdy, Smo, Syr, Shr, Smn
   caldat, EndJul, Edy, Emo, Eyr, Ehr, Emn
   
   Xrng = [StartJul, EndJul]
      
   Fn = out_fn_core + string(Shr, Smn, Ehr, Emn, format = '("_",i02,i02,"_",i02,i02)')
   
   WhFlux1 = where(Flux1.JdCor ne 0.0)
   if keyword_set(FLux2) then WhFlux2 = where(Flux2.JdCor ne 0.0)
   WhMag = where(magJday ne 0.0)   

   FluxColor = [color_scale[2], color_scale[5], color_scale[8], color_scale[9], $
      color_scale[6], color_scale[12], color_scale[10], color_scale[7], color_scale[0]]
   PALabScale = [45,85,105,65,125,25,165,145,5] ; deg (inverted)

   plot, Flux1.JdCor[0,WhFlux1], Flux1.Flux[0,WhFlux1], /nodata, $ 
      xr = Xrng, xstyle=1, xtickformat='xticksMAGED', xticks = nxticks, xminor = nxminor, $
      ytitle = 'cm!U-2!Ns!U-1!Nsr!U-1!NkeV!U-1!N', /ylog, yr = Yrng, ystyle = 1, $
      charsize = ChSz, title = DateStr + ' ' + flux1_title

   for i = 0, 8 do begin
      oplot, Flux1.JdCor[i,WhFlux1], Flux1.Flux[i,WhFlux1], color = FluxColor[i]
;      xyouts, EndJul + 0.005*(EndJul-StartJul), PALabScale[i], $
;         string(i+1, format = '(i1)'), color = FluxColor[i], charsize = 0.9
   endfor

   if keyword_set(Flux2) then begin
      plot, Flux2.JdCor[0,WhFlux2], Flux2.Flux[0,WhFlux2], /nodata, $ 
         xr = Xrng, xstyle=1, xtickformat='xticksMAGED', xticks = nxticks, xminor = nxminor, $
         ytitle = 'cm!U-2!Ns!U-1!Nsr!U-1!NkeV!U-1!N', /ylog, yr = Yrng2, ystyle = 1, $
         charsize = ChSz, title = flux2_title

      for i = 0, 8 do begin
         oplot, Flux2.JdCor[i,WhFlux2], Flux2.Flux[i,WhFlux2], color = FluxColor[i]
;      xyouts, EndJul + 0.005*(EndJul-StartJul), PALabScale[i], $
;         string(i+1, format = '(i1)'), color = FluxColor[i], charsize = 0.9
      endfor
   endif
   
   plot, magJday[WhMag], magArray[WhMag,0], /nodata, $
      xr = Xrng, xstyle=1, xtickformat='xticksMAGED', xticks = nxticks, xminor = nxminor, $
      ytitle = 'B (nT)', yr = [-20, 120], ystyle = 1, yticks = 7, yminor = 5, $
      charsize = ChSz, title = mag_title

   HstrArr = ["Hp","He","Hn","Ht"]
   HlabScale = [60, 30, 0, 90]
   for i = 0, 3 do begin
      oplot, magJday[WhMag], magArray[WhMag,i], color = FluxColor[i]
      xyouts, EndJul + 0.005*(EndJul-StartJul), HLabScale[i], $
         HstrArr[i], color = FluxColor[i], charsize = 0.9
   endfor
   
   if plot_type eq 'png' then begin
   Img = TVRD(0, TRUE = 1)
   PlotnamePNG = Fn + '.png'
   WRITE_PNG, PlotnamePNG, Img
   endif else if plot_type eq 'ps' then device, /close

end



;----------------------------------------------------------------------------------

pro plot_fit_eflux, N1min, MAGEDEnergy, GEfactor, jflux, MAGEDRMQ, MAGEDRMCovQQ, $
   MAGEDPLQ, MAGEDPLCovQQ, MAGEDEEQ, MAGEDEECovQQ, MAGEDDMQ, MAGEDDMCovQQ, wk, origy, $
   electronFluxes, eLambda, elambda_pl, elambda_rm, elambda_dm, elambda_ee, eSigma

; Modified Dec 3, 2010

   common weighting, dt, energy, Kwfe 
   
   mass_e_keV = 510.998910d0   
   plot_energy = energy
   energy_coeff = plot_energy*(plot_energy + 2.0d0*mass_e_keV)/2.0d0/mass_e_keV
   electronFluxes = []
   eLambda = []
   eLambda_pl = []
   elambda_ee = []
   elambda_rm = []
   elambda_dm = []
   eSigma = []
   
   for i = 0, N1min-1 do begin

     if MAGEDRMQ[i,0] ne -1.0e31 then begin ;this was in plot_fit_pflux
         covqqRM = reform(MAGEDRMCovQQ[i,*,*])
         QsRM = reform(MAGEDRMQ[i,*])
         RM = energy_coeff*exp(QsRM[0] + QsRM[1]*plot_energy)
         logflux_standard_deviation_2, QsRM, covqqRM, plot_energy, $
            'relmaxwell', 'relmaxwell_gradient', logRMsd, RMsd

         covqqPL = reform(MAGEDPLCovQQ[i,*,*])
         QsPL = reform(MAGEDPLQ[i,*])
         PL = exp(QsPL[0] - QsPL[1]*alog(plot_energy))
         logflux_standard_deviation_2, QsPL, covqqPL, plot_energy, $
            'powerlaw', 'powerlaw_gradient', logPLsd, PLsd

         covqqEE = reform(MAGEDEECovQQ[i,*,*])
         QsEE = reform(MAGEDEEQ[i,*])
         EE = exp(QsEE[0] + QsEE[1]*plot_energy)
         logflux_standard_deviation_2, QsEE, covqqEE, plot_energy, $
            'energyexponential', 'energyexponential_gradient', logEEsd, EEsd

         covqqDM = reform(MAGEDDMCovQQ[i,*,*])
         QsDM = reform(MAGEDDMQ[i,*])
         DM = energy_coeff*(exp(QsDM[0] + QsDM[1]*plot_energy)+exp(QsDM[2] + QsDM[3]*plot_energy))
         logflux_standard_deviation_4, QsDM, covqqDM, plot_energy, $
            'fdrelmaxwell', 'fdrelmaxwell_gradient', logDMsd, DMsd

         ln_combined_flux = (wk[i,0]*alog(RM) + wk[i,1]*alog(PL) + wk[i,2]*alog(EE) + wk[i,3]*alog(DM))
         combined_flux = exp(ln_combined_flux)      
         
         ln_rm_flux = alog(rm)
         rm_flux = exp(ln_rm_flux)
         
         ln_pl_flux = alog(PL)
         pl_flux = exp(ln_pl_flux)
         
         ln_ee_flux = alog(ee)
         ee_flux = exp(ln_ee_flux)

         ln_dm_flux = alog(dm)
         dm_flux = exp(ln_dm_flux)
         
	;if i eq 2790 then stop
         sigma_combined_logflux_obrien = sqrt(wk[i,0]*(logRMsd^2. + alog(RM)^2.) + $
                                              wk[i,1]*(logPLsd^2. + alog(PL)^2.) + $
                                              wk[i,2]*(logEEsd^2. + alog(EE)^2.) + $
					      wk[i,3]*(logDMsd^2. + alog(DM)^2.) - ln_combined_flux^2.) 
         sigma_combined_flux_obrien = combined_flux * sigma_combined_logflux_obrien
         sigma_combined_logflux_me =     sqrt((wk[i,0]*logRMsd)^2. + $
                                              (wk[i,1]*logPLsd)^2. + $
                                              (wk[i,2]*logEEsd)^2. + $
					      (wk[i,3]*logDMsd)^2.) 
         sigma_combined_flux_me = combined_flux * sigma_combined_logflux_me
	 sigma_combined_logflux_ethan =  sqrt((wk[i,0]*logRMsd^2.) + $
					      (wk[i,1]*logPLsd^2.) + $
					      (wk[i,2]*logEEsd^2.) + $
					      (wk[i,3]*logDMsd^2.))
	 sigma_combined_flux_ethan = combined_flux * sigma_combined_logflux_ethan
              
              
                  
	 lambda = simpleForwardModel(Kwfe,combined_flux,dt,energy)
	 
	 lambda_pl = lambda
	 lambda_ee = lambda
	 lambda_rm = lambda
	 lambda_dm = lambda
	 
	 lambda_pl = simpleForwardModel(kwfe, pl_flux, dt, energy)
	 lambda_ee = simpleForwardModel(kwfe, ee_flux, dt, energy)
	 lambda_rm = simpleForwardModel(kwfe, rm_flux, dt, energy)
	 lambda_dm = simpleForwardModel(kwfe, dm_flux, dt, energy)
	 
     endif else begin
         combined_flux = fltarr(N_ELEMENTS(energy))-1.e-31
         lambda = fltarr(4)-1.e-31
         lamda_pl = fltarr(4)-1.e-31
	 sigma_combined_flux_ethan = fltarr(N_ELEMENTS(energy))-1.e-31
	 sigma_combined_flux_me = fltarr(N_ELEMENTS(energy))-1.e-31
	 sigma_combined_flux_obrien = fltarr(N_ELEMENTS(energy))-1.e-31
     endelse
       electronFluxes = [[electronFluxes],[combined_flux]]
       eLambda = [[eLambda],[lambda]]
       
       ; Inserted to deal with stupid undefined issues occuring for Paper 5 (investigate further)
       
       lambda_pl = lambda
       lambda_ee = lambda
       lambda_rm = lambda
       lambda_dm = lambda
       
       elambda_pl = [[elambda_pl], [lambda_pl]]
       elambda_ee = [[elambda_ee], [lambda_ee]]
       elambda_rm = [[elambda_rm], [lambda_rm]]
       elambda_dm = [[elambda_dm], [lambda_dm]]
       
       ;eSigma = [[eSigma],[sigma_combined_flux_ethan]]
       eSigma = [[eSigma],[sigma_combined_flux_obrien]]
   endfor
; stop
end

;----------------------------------------------------------------------------------

pro plot_fit_pflux, N1min, MAGPDEnergy, GEfactor, jflux, jflux_energy, MAGPDRMQ, MAGPDRMCovQQ, $
   MAGPDPLQ, MAGPDPLCovQQ, MAGPDEEQ, MAGPDEECovQQ, MAGPDDMQ, MAGPDDMCovQQ, wkp, origy,protonFluxes,pLambda,pSigma

   common weighting, dt, energy, Kwfe, Kwfp
   
   mass_p_keV = 938272.013d0 
   plot_energy = energy
   energy_coeff = plot_energy*(plot_energy + 2.0d0*mass_p_keV)/2.0d0/mass_p_keV
   protonFluxes = []
   pLambda = []
   pSigma = []
   
   for i = 0, N1min-1 do begin
         if MAGPDRMQ[i,0] ne -1.0e31 then begin
  
         covqqRM = reform(MAGPDRMCovQQ[i,*,*])
         QsRM = reform(MAGPDRMQ[i,*])
         RM = energy_coeff*exp(QsRM[0] + QsRM[1]*plot_energy)
         logflux_standard_deviation_2, QsRM, covqqRM, plot_energy, $
            'relmaxwellp', 'relmaxwellp_gradient', logRMsd, RMsd

         covqqPL = reform(MAGPDPLCovQQ[i,*,*])
         QsPL = reform(MAGPDPLQ[i,*])
         PL = exp(QsPL[0] - QsPL[1]*alog(plot_energy))
         logflux_standard_deviation_2, QsPL, covqqPL, plot_energy, $
            'powerlaw', 'powerlaw_gradient', logPLsd, PLsd

         covqqEE = reform(MAGPDEECovQQ[i,*,*])
         QsEE = reform(MAGPDEEQ[i,*])
         EE = exp(QsEE[0] + QsEE[1]*plot_energy)
         logflux_standard_deviation_2, QsEE, covqqEE, plot_energy, $
            'energyexponential', 'energyexponential_gradient', logEEsd, EEsd

         covqqDM = reform(MAGPDDMCovQQ[i,*,*])
         QsDM = reform(MAGPDDMQ[i,*])
         DM = energy_coeff*(exp(QsDM[0] + QsDM[1]*plot_energy)+exp(QsDM[2] + QsDM[3]*plot_energy))
         logflux_standard_deviation_4, QsDM, covqqDM, plot_energy, $
            'fdrelmaxwellp', 'fdrelmaxwellp_gradient', logDMsd, DMsd
 
         ln_combined_flux = (wkp[i,0]*alog(RM) + wkp[i,1]*alog(PL) + wkp[i,2]*alog(EE) + wkp[i,3]*alog(DM))
         combined_flux = exp(ln_combined_flux)
         sigma_combined_logflux_obrien = sqrt(wkp[i,0]*(logRMsd^2. + alog(RM)^2.) + $
                                              wkp[i,1]*(logPLsd^2. + alog(PL)^2.) + $
                                              wkp[i,2]*(logEEsd^2. + alog(EE)^2.) + $
					      wkp[i,3]*(logDMsd^2. + alog(DM)^2.) - ln_combined_flux^2.) 
         sigma_combined_flux_obrien = combined_flux * sigma_combined_logflux_obrien
         sigma_combined_logflux_me =     sqrt((wkp[i,0]*logRMsd)^2. + $
                                              (wkp[i,1]*logPLsd)^2. + $
                                              (wkp[i,2]*logEEsd)^2. + $
					      (wkp[i,3]*logDMsd)^2.) 
         sigma_combined_flux_me = combined_flux * sigma_combined_logflux_me
	 sigma_combined_logflux_ethan =  sqrt((wkp[i,0]*logRMsd^2.) + $
					      (wkp[i,1]*logPLsd^2.) + $
					      (wkp[i,2]*logEEsd^2.) + $
					      (wkp[i,3]*logDMsd^2.))
	 sigma_combined_flux_ethan = combined_flux * sigma_combined_logflux_ethan
                  

	lambda = simpleForwardModel(Kwfp,combined_flux,dt,energy)
         endif else begin
	   combined_flux = fltarr(N_ELEMENTS(energy))-1.e-31
           lambda = fltarr(6)-1.e-31
	   sigma_combined_flux_ethan = fltarr(N_ELEMENTS(energy))-1.e-31
	   sigma_combined_flux_obrien = fltarr(N_ELEMENTS(energy))-1.e-31
	   sigma_combined_flux_me = fltarr(N_ELEMENTS(energy))-1.e-31
         endelse
      protonFluxes = [[protonFluxes],[combined_flux]]
      pLambda = [[pLambda],[lambda]]
      ;pSigma = [[pSigma],[sigma_combined_flux_ethan]]
      pSigma = [[pSigma],[sigma_combined_flux_obrien]]

   endfor
   
end

;----------------------------------------------------------------------------------

pro outputNCfile,energy,eflux,pflux,eSigma,pSigma,eLambda,eLambda_pl,eLambda_dm,eLambda_rm,eLambda_ee,pLambda,geogLL,foflLL,time,MLT,lValue,pitch,electronCount,protonCount,new_electronCount,eRMQ,ePLQ,eEEQ,eDMQ,pRMQ,pPLQ,pEEQ,pDMQ,w,wp,sat,satName,date,dt,angle,msec_time,protoncontamination,bfofl,blocal,blc_angle
;this procedure will create a NetCDF output of the corrected MEPED flux data. 
;sat is a string with the short name of the satellite being analyzed, so 'm02' is metop2.
;satName is a string the the full name of the satellite being analyzed, so 'METOP02' is metop2.
;			stop
			;calculate real time as opposed to epoch time, because I am human and do not think in epoch time.
			rtime = dblarr(N_ELEMENTS(time))
			;for i = 0, N_ELEMENTS(time)-1 do begin
	    ;                   CDF_EPOCH, time[i], yr, mo, dy, hr, mn, sc, milli, /BREAK
      ;  	               rTime[i] = ((milli/1000.+sc)/60.+mn)/60.+hr
			;endfor
			;stop
                        ;output netcdf file for day of measurements
                        Ecounts = eflux
                        EcountsATTR = 'Corrected Electron Flux Rate'
                        EcountsUnit = 'N/cm^2/sr/keV'

                        Pcounts = pflux
                        PcountsATTR = 'Corrected Proton Flux Rate'
                        PcountsUnit = 'N/cm^2/sr/keV'

                        Eerror = eSigma
                        EerrorATTR = 'Error on Corrected Electron Count Rate, 1-sigma'
                        EerrorUnit = 'N/cm^2/sr/keV'

                        Perror = pSigma
                        PerrorATTR = 'Error on Corrected Proton Count Rate, 1-sigma'
                        PerrorUnit = 'N/cm^2/sr/keV'

                        geog_lat = REFORM(geogLL[0,*])
                        geog_latATTR = 'Satellite Geographic Latitude'
                        geog_latUnit = 'degrees -90 to 90'

                        geog_lon = REFORM(geogLL[1,*])
                        geog_lonATTR = 'Satellite Geographic Longitude'
                        geog_lonUnit = 'degrees 0 to 360'

                        fofl_lat = REFORM(foflLL[0,*])
                        fofl_latATTR = 'Satellite Foot of the Field Line Geographic Latitude'
                        fofl_latUnit = 'degrees -90 to 90'

                        fofl_lon = REFORM(foflLL[1,*])
                        fofl_lonATTR = 'Satellite Foot of the Field Line Geographic Longitude'
                        fofl_lonUnit = 'degrees 0 to 360'

                        ;convert time from epoch to something useful here.
                        timeATTR = 'Time of measurement'
                        timeUnit = 'Hours'

                 ;       if rtime[0] gt 23. then rtime[0] = rtime[0] - 24. ;check for yesterday leaking into today
                 ;       rtimecheck = where((rtime[1:-1] - rtime[0:-2]) lt 0.,/null) ;check for today leaking into tomorrow
                 ;       if rtimecheck ne !NULL then rtime[rtimecheck+1] = rtime[rtimecheck+1]+24.

                 ;       rTimeATTR = 'Time of measurment'
                 ;       rTimeUnit = 'Seconds'

                        MLTATTR = 'Satellite Magnetic Local Time'
                        MLTUnit = 'degrees 0 to 360'

                        lValueATTR = 'Satellite L-Value'
                        lValueUnit = ' '
			
						pitch = REFORM(pitch)
                        pitchATTR = 'Pitch Angle'
                        pitchUnit = 'degrees -180 to 180'
                        
                        Bfoflattr = 'Magnetic Field at measurement'
                        Bfoflunit = 'Tesla'

                        Blocalattr = 'Magnetic field at satellite'
                        blocalunit = 'Tesla'

                        BLC_angleattr = 'BLC Angle at Satellite'
                        BLC_angleunit = 'degrees -180 to 180'

                        energyATTR = 'Center Energy of Electron/Proton Fluxes'
                        energyUnit = 'KeV'

                        ;edgesATTR = 'Edge Energies of Electron/Proton Count Rate'
                        ;edgesDescription = 'Fields 1 and 2 coordinate to energy 1. Fields 2 and 3 coordinate to energy 2. etc.'
                        ;edgesUnit = 'KeV'

                        ;Emask = mask[0:26,*]
                        ;EmaskATTR = 'Electron Count Rate Mask'
                        ;EmaskDescription = 'Created for count rates where error is equal to or greater than 100% of the count rate'
                        ;EmaskUnit = 'KeV'

                        ;Pmask = mask[27:-1,*]
                        ;PmaskATTR = 'Proton Count Rate Mask'
                        ;PmaskDescription = 'Created for count rates where error is equal to or greater than 100% of the count rate'
                        ;PmaskUnit = 'KeV'

                        EOcounts = electronCount
                        EOcountsATTR = 'Original Electron Count Rate'
                        EOcountsUnit = 'counts/sec'

                        POcounts = protonCount
                        POcountsATTR = 'Original Proton Count Rate'
                        POcountsUnit = 'counts/sec'

                        EOcounts_corrected = new_electronCount
                        EOcounts_correctedATTR = 'Corrected Electron Count Rate'
                        EOcounts_correctedUnit = 'counts/sec'

                        EOcounts_Lambda = eLambda/dt
                        EOcounts_LambdaATTR = 'Forward Model Corrected Electron Count Rate'
                        EOcounts_LambdaUnit = 'counts/sec'
                        
                        EOcounts_Lambda_pl = eLambda_pl/dt
                        EOcounts_LambdaATTR_pl = 'Forward Model Corrected Electron Count Rate for Powerlaw'
                        EOcounts_LambdaUnit_pl = 'counts/sec'

                        EOcounts_Lambda_dm = eLambda_dm/dt
                        EOcounts_LambdaATTR_dm = 'Forward Model Corrected Electron Count Rate for Double Maxwellian'
                        EOcounts_LambdaUnit_dm = 'counts/sec'

                        EOcounts_Lambda_rm = eLambda_rm/dt
                        EOcounts_LambdaATTR_rm = 'Forward Model Corrected Electron Count Rate for Relativistic Maxwellian'
                        EOcounts_LambdaUnit_rm = 'counts/sec'

                        EOcounts_Lambda_ee = eLambda_ee/dt
                        EOcounts_LambdaATTR_ee = 'Forward Model Corrected Electron Count Rate for Energy Exponential'
                        EOcounts_LambdaUnit_ee = 'counts/sec'

                        POcounts_Lambda = pLambda/dt
                        POcounts_LambdaATTR = 'Forward Model Proton Count Rate'
                        POcounts_LambdaUnit = 'counts/sec'
                        
                        Proton_Contamination = ProtonContamination
                        Proton_ContaminationATTR = 'Proton Contamination'
                        Proton_ContaminationUnit = 'counts/sec'

			eRMq = TRANSPOSE(eRMq)
			eRMqATTR1 = 'Parameters for Electron Relativistic Maxwellian Spectrum'
			eRMqATTR2 = 'f(E) = E*(1 + E/E0/2)*exp(q1 + q2*E)'
			eRMqATTR3 = 'E0 = 511 keV'

			pRMq = TRANSPOSE(pRMq)
			pRMqATTR1 = 'Parameters for Proton Relativistic Maxwellian Spectrum'
			pRMqATTR2 = 'f(E) = E*(1 + E/E0/2)*exp(q1 + q2*E)'
			pRMqATTR3 = 'E0 = 938 MeV'

			ePLq = TRANSPOSE(ePLq)
			ePLqATTR1 = 'Parameters for Electron Power Law Spectrum'
			ePLqATTR2 = 'f(E) = exp(q1 - q2*ln(E))'

			pPLq = TRANSPOSE(pPLq)
			pPLqATTR1 = 'Parameters for Proton Power Law Spectrum'
			pPLqATTR2 = 'f(E) = exp(q1 - q2*ln(E))'

			eEEq = TRANSPOSE(eEEq)
			eEEqATTR1 = 'Parameters for Electron Energy Exponential Spectrum'
			eEEqATTR2 = 'f(E) = exp(q1 + q2*E)'

			pEEq = TRANSPOSE(pEEq)
			pEEqATTR1 = 'Parameters for Proton Energy Exponential Spectrum'
			pEEqATTR2 = 'f(E) = exp(q1 + q2*E)'

			eDMq = TRANSPOSE(eDMq)
			eDMqATTR1 = 'Parameters for Electron Double Relativistic Maxwellian Spectrum'
			eDMqATTR2 = 'f(E) = E*(1 + E/E0/2)*[exp(q1 + q2*E)+exp(q3 + q4*E)]'
			eDMqATTR3 = 'E0 = 511 keV'

			pDMq = TRANSPOSE(pDMq)
			pDMqATTR1 = 'Parameters for Proton Double Relativistic Maxwellian Spectrum'
			pDMqATTR2 = 'f(E) = E*(1 + E/E0/2)*[exp(q1 + q2*E)+exp(q3 + q4*E)]'
			pDMqATTR3 = 'E0 = 938 MeV'

			w = TRANSPOSE(w)
			wATTR1 = 'Weighting of Log Electron Fluxes'
			wATTR2 = 'Relativistic Maxwellian, Power Law, Exponential, and Double Relativistic Maxwellian'
			wATTR3 = 'Multiply respective elements to ln of spectra, sum products together, and exp to get combined spectrum.'

			wp = TRANSPOSE(wp)
			wpATTR1 = 'Weighting of Log Proton Fluxes'
			wpATTR2 = 'Relativistic Maxwellian, Power Law, Exponential, and Double Relativistic Maxwellian'
			wpATTR3 = 'Multiply respective elements to ln of spectra, sum products together, and exp to get combined spectrum.'


      ;stop
                        dataType = string(dt,format='(I2.2)')+' second average data'
                        ofilename = 'POES_combinedSpectrum_'+sat+'_'+string(angle,FORMAT='(I2.2)')+'_'+date[0]+date[1]+date[2]+'.nc'
                        ;stop   
                        ;ofile = NCDF_CREATE('Z:\MEE_Data\Processed_MEPED_Data\Kathy\'+ofilename,/CLOBBER)
                        ;ofile = NCDF_CREATE('Z:\\MEE_Data\Processed_MEPED_Data\NOAA15\'+ofilename,/CLOBBER)
                        ofile = NCDF_CREATE('D:\POES_Data\MPE_Software\Processed_MEPED_Data\'+satName+'\'+ofilename,/CLOBBER)
                        ;ofile = NCDF_CREATE('Z:\MEE_Data\Allison_Data\'+ofilename, /CLOBBER)
                        ;define dimensions
                        energydimID = NCDF_DIMDEF(ofile,'energy',N_ELEMENTS(energy))
                        ;edgesdimID = NCDF_DIMDEF(ofile,'edges',N_ELEMENTS(edges))
                        PtelescopedimID = NCDF_DIMDEF(ofile,'proton_telescopes',6)
                        EtelescopedimID = NCDF_DIMDEF(ofile,'electron_telescopes',3)
                        Etelescope_correcteddimID = NCDF_DIMDEF(ofile,'electron_telescopes_and_E4',4)
                  			parameters2ID = NCDF_DIMDEF(ofile,'2_Parameter_Spectrum',2)
                  			parameters4ID = NCDF_DIMDEF(ofile,'4_Parameter_Spectrum',4)
                  			spectraID = NCDF_DIMDEF(ofile,'number_of_spectra',4)
                        timedimID = NCDF_DIMDEF(ofile,'time',/UNLIMITED)
                        ;define variables
                       ; stop
                        timeID = NCDF_VARDEF(ofile,'time',timedimID,/DOUBLE)
                        ;msectimeid = ncdf_vardef(ofile, 'msec_time', timedimID, /DOUBLE)
                       ;rtimeID = NCDF_VARDEF(ofile,'rtime',timedimID,/DOUBLE)
                        geoglatID = NCDF_VARDEF(ofile,'geogLat',timedimID,/DOUBLE)
                        geoglonID = NCDF_VARDEF(ofile,'geogLon',timedimID,/DOUBLE)
                        fofllatID = NCDF_VARDEF(ofile,'foflLat',timedimID,/DOUBLE)
                        fofllonID = NCDF_VARDEF(ofile,'foflLon',timedimID,/DOUBLE)
                        MLTID = NCDF_VARDEF(ofile,'MLT',timedimID,/DOUBLE)
                        lvalueID = NCDF_VARDEF(ofile,'lValue',timedimID,/DOUBLE)
                        energyID = NCDF_VARDEF(ofile,'energy',energydimID,/FLOAT)
                        ;edgesID = NCDF_VARDEF(ofile,'edges',edgesdimID,/FLOAT)
                        EOcountsID = NCDF_VARDEF(ofile,'EOcounts',[etelescopedimID,timedimID],/DOUBLE)
                        ProtonContaminationID = NCDF_VARDEF(ofile,'Proton_Contamination',[etelescopedimID,timedimID],/DOUBLE)
                        EOcounts_LambdaID = NCDF_VARDEF(ofile,'EOcounts_Lambda',[etelescope_CORRECTEDdimID,timedimID],/DOUBLE)
                        EOcounts_LambdaID_pl = NCDF_VARDEF(ofile,'EOcounts_Lambda_pl',[etelescope_CORRECTEDdimID,timedimID],/DOUBLE)
                        EOcounts_LambdaID_dm = NCDF_VARDEF(ofile,'EOcounts_Lambda_dm',[etelescope_CORRECTEDdimID,timedimID],/DOUBLE)
                        EOcounts_LambdaID_rm = NCDF_VARDEF(ofile,'EOcounts_Lambda_rm',[etelescope_CORRECTEDdimID,timedimID],/DOUBLE)
                        EOcounts_LambdaID_ee = NCDF_VARDEF(ofile,'EOcounts_Lambda_ee',[etelescope_CORRECTEDdimID,timedimID],/DOUBLE)
                        EOcounts_correctedID = NCDF_VARDEF(ofile,'EOcounts_corrected',[etelescope_correcteddimID,timedimID],/DOUBLE)
                        POcountsID = NCDF_VARDEF(ofile,'POcounts',[ptelescopedimID,timedimID],/DOUBLE)
                        POcounts_LambdaID = NCDF_VARDEF(ofile,'POcounts_Lambda',[ptelescopedimID,timedimID],/DOUBLE)
                        EcountsID = NCDF_VARDEF(ofile,'Ecounts',[energydimID,timedimID],/DOUBLE)
                        PcountsID = NCDF_VARDEF(ofile,'Pcounts',[energydimID,timedimID],/DOUBLE)
                        EerrorID = NCDF_VARDEF(ofile,'Eerror',[energydimID,timedimID],/DOUBLE)
                        PerrorID = NCDF_VARDEF(ofile,'Perror',[energydimID,timedimID],/DOUBLE)
                        ;EmaskID = NCDF_VARDEF(ofile,'Emask',[energydimID,timedimID],/DOUBLE)
                        ;PmaskID = NCDF_VARDEF(ofile,'Pmask',[energydimID,timedimID],/DOUBLE)
                        ;EprioriID = NCDF_VARDEF(ofile,'Epriori',[energydimID,timedimID],/DOUBLE)
                        ;PprioriID = NCDF_VARDEF(ofile,'Ppriori',[energydimID,timedimID],/DOUBLE)
                  			ERMqID = NCDF_VARDEF(ofile,'ERMq',[parameters2ID,timedimID],/DOUBLE)
                  			EPLqID = NCDF_VARDEF(ofile,'EPLq',[parameters2ID,timedimID],/DOUBLE)
                  			EEEqID = NCDF_VARDEF(ofile,'EEEq',[parameters2ID,timedimID],/DOUBLE)
                  			EDMqID = NCDF_VARDEF(ofile,'EDMq',[parameters4ID,timedimID],/DOUBLE)
                  			PRMqID = NCDF_VARDEF(ofile,'PRMq',[parameters2ID,timedimID],/DOUBLE)
                  			PPLqID = NCDF_VARDEF(ofile,'PPLq',[parameters2ID,timedimID],/DOUBLE)
                  			PEEqID = NCDF_VARDEF(ofile,'PEEq',[parameters2ID,timedimID],/DOUBLE)
                  			PDMqID = NCDF_VARDEF(ofile,'PDMq',[parameters4ID,timedimID],/DOUBLE)
                  			WID = NCDF_VARDEF(ofile,'W',[spectraID,timedimID],/DOUBLE)
                  			WpID = NCDF_VARDEF(ofile,'Wp',[spectraID,timedimID],/DOUBLE)
                        pitchID = NCDF_VARDEF(ofile,'pitch',timedimID,/DOUBLE)
                        BLC_FootID = NCDF_VARDEF(ofile,'Bfofl',timedimID,/DOUBLE)
                        BLC_SatID = NCDF_VARDEF(ofile,'Blocal',timedimID,/DOUBLE)
                        BLC_AngleID = NCDF_VARDEF(ofile,'BLC_Angle',timedimID,/DOUBLE)

                        ;define attributes
                        NCDF_ATTPUT,ofile,timeID,'Full_Name',timeATTR,/CHAR
                        NCDF_ATTPUT,ofile,timeID,'Units',timeUnit,/CHAR
                      ; NCDF_ATTPUT,ofile,rtimeID,'Full_Name',rtimeATTR,/CHAR
                      ; NCDF_ATTPUT,ofile,rtimeID,'Units',rtimeUnit,/CHAR
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
                        NCDF_ATTPUT,ofile,EOcountsID,'Full_Name',EOcountsATTR,/CHAR
                        NCDF_ATTPUT,ofile,EOcountsID,'Units',EOcountsUnit,/CHAR
                        NCDF_ATTPUT,ofile,EOcounts_LambdaID,'Full_Name',EOcounts_LambdaATTR,/CHAR
                        NCDF_ATTPUT,ofile,EOcounts_LambdaID,'Units',EOcounts_LambdaUnit,/CHAR
                        NCDF_ATTPUT,ofile,EOcounts_correctedID,'Full_Name',EOcounts_correctedATTR,/CHAR
                        NCDF_ATTPUT,ofile,EOcounts_correctedID,'Units',EOcounts_correctedUnit,/CHAR
                        NCDF_ATTPUT,ofile,POcountsID,'Full_Name',POcountsATTR,/CHAR
                        NCDF_ATTPUT,ofile,POcountsID,'Units',POcountsUnit,/CHAR
                        NCDF_ATTPUT,ofile,ProtonContaminationID,"Bad_Data_Value",-999.,/DOUBLE
                        NCDF_ATTPUT,ofile,POcounts_LambdaID,'Full_Name',POcounts_LambdaATTR,/CHAR
                        NCDF_ATTPUT,ofile,POcounts_LambdaID,'Units',POcounts_LambdaUnit,/CHAR
                        NCDF_ATTPUT,ofile,EcountsID,'Full_Name',EcountsATTR,/CHAR
                        NCDF_ATTPUT,ofile,EcountsID,'Units',EcountsUnit,/CHAR
                        NCDF_ATTPUT,ofile,EcountsID,'Bad_Data_Value',-999.,/DOUBLE
                        NCDF_ATTPUT,ofile,PcountsID,'Full_Name',PcountsATTR,/CHAR
                        NCDF_ATTPUT,ofile,PcountsID,'Units',PcountsUnit,/CHAR
                        NCDF_ATTPUT,ofile,PcountsID,'Bad_Data_Value',-999.,/DOUBLE
                        NCDF_ATTPUT,ofile,EerrorID,'Full_Name',EerrorATTR,/CHAR
                        NCDF_ATTPUT,ofile,EerrorID,'Units',EerrorUnit,/CHAR
                        NCDF_ATTPUT,ofile,PerrorID,'Full_Name',PerrorATTR,/CHAR
                        NCDF_ATTPUT,ofile,PerrorID,'Units',PerrorUnit,/CHAR
            						NCDF_ATTPUT,ofile,ERMqID,'Full_Name',eRMqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,ERMqID,'Spectrum_Equation',eRMqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,ERMqID,'Rest_Energy',eRMqATTR3,/CHAR
            						NCDF_ATTPUT,ofile,EPLqID,'Full_Name',ePLqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,EPLqID,'Spectrum_Equation',ePLqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,EEEqID,'Full_Name',eEEqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,EEEqID,'Spectrum_Equation',eEEqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,EDMqID,'Full_Name',eDMqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,EDMqID,'Spectrum_Equation',eDMqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,EDMqID,'Rest_Energy',eDMqATTR3,/CHAR
            						NCDF_ATTPUT,ofile,PRMqID,'Full_Name',pRMqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,PRMqID,'Spectrum_Equation',pRMqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,PRMqID,'Rest_Energy',pRMqATTR3,/CHAR
            						NCDF_ATTPUT,ofile,PPLqID,'Full_Name',pPLqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,PPLqID,'Spectrum_Equation',pPLqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,PEEqID,'Full_Name',pEEqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,PEEqID,'Spectrum_Equation',pEEqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,PDMqID,'Full_Name',pDMqATTR1,/CHAR
            						NCDF_ATTPUT,ofile,PDMqID,'Spectrum_Equation',pDMqATTR2,/CHAR
            						NCDF_ATTPUT,ofile,PDMqID,'Rest_Energy',pDMqATTR3,/CHAR
            						NCDF_ATTPUT,ofile,WID,'Full_Name',wATTR1,/CHAR
            						NCDF_ATTPUT,ofile,WID,'Order_of_elements',wATTR2,/CHAR
            						NCDF_ATTPUT,ofile,WID,'How_to_use',wATTR3,/CHAR
            						NCDF_ATTPUT,ofile,WpID,'Full_Name',wpATTR1,/CHAR
            						NCDF_ATTPUT,ofile,WpID,'Order_of_elements',wpATTR2,/CHAR
            						NCDF_ATTPUT,ofile,WpID,'How_to_use',wpATTR3,/CHAR
                        NCDF_ATTPUT,ofile,pitchID,'Full_Name',pitchATTR,/CHAR
                        NCDF_ATTPUT,ofile,pitchID,'Units',pitchUnit,/CHAR
                        NCDF_ATTPUT,ofile,BLC_FootID,'Full_Name',BLC_FootID,/CHAR
                        NCDF_ATTPUT,ofile,BLC_FootID,'Units',BLC_FootID,/CHAR
                        NCDF_ATTPUT,ofile,BLC_SatID,'Full_Name',BLC_SatID,/CHAR
                        NCDF_ATTPUT,ofile,BLC_SatID,'Units',BLC_SatID,/CHAR
                        NCDF_ATTPUT,ofile,BLC_angleID,'Full_Name',BLC_angleID,/CHAR
                        NCDF_ATTPUT,ofile,BLC_angleID,'Units',BLC_angleID,/CHAR
                        NCDF_ATTPUT,ofile,'Data_Type',dataType,/CHAR,/GLOBAL
                        NCDF_ATTPUT,ofile,'Author','Joshua Pettit, University of Colorado Boulder, LASP',/CHAR,/GLOBAL
                        NCDF_ATTPUT,ofile,'Date_Created','Spring 2020',/CHAR,/GLOBAL
                        NCDF_ATTPUT,ofile,'Satellite',satName,/CHAR,/GLOBAL

                        ;input variables
                        NCDF_CONTROL,ofile,/ENDEF
                      ;  stop
                       ; time = (time/1000.)-(time[0]/1000.)
                      ;  rtime = rtime-rtime[0]
                        NCDF_VARPUT,ofile,timeID,time
                      ;  NCDF_VARPUT,ofile,rtimeID,rtime
                     ;   NCDF_VARPUT,ofile,msectimeID,msec_time
                        NCDF_VARPUT,ofile,geoglatID,geog_lat
                        NCDF_VARPUT,ofile,geoglonID,geog_lon
                        NCDF_VARPUT,ofile,fofllatID,fofl_lat
                        NCDF_VARPUT,ofile,fofllonID,fofl_lon
                        NCDF_VARPUT,ofile,MLTID,MLT
                        NCDF_VARPUT,ofile,lvalueID,lvalue
                        NCDF_VARPUT,ofile,energyID,energy
                        NCDF_VARPUT,ofile,EOcountsID,EOcounts
                        NCDF_VARPUT,ofile,EOcounts_LambdaID,EOcounts_Lambda
                        NCDF_VARPUT,ofile,EOcounts_LambdaID_pl,EOcounts_Lambda_pl
                        NCDF_VARPUT,ofile,EOcounts_LambdaID_dm,EOcounts_Lambda_dm
                        NCDF_VARPUT,ofile,EOcounts_LambdaID_rm,EOcounts_Lambda_rm
                        NCDF_VARPUT,ofile,EOcounts_LambdaID_ee,EOcounts_Lambda_ee
                        NCDF_VARPUT,ofile,EOcounts_correctedID,EOcounts_corrected
                        NCDF_VARPUT,ofile,POcountsID,POcounts
                        NCDF_VARPUT,ofile,POcounts_LambdaID,POcounts_Lambda
                        NCDF_VARPUT,ofile,ProtonContaminationID,Proton_Contamination
                        NCDF_VARPUT,ofile,EcountsID,Ecounts
                        NCDF_VARPUT,ofile,PcountsID,Pcounts
                        NCDF_VARPUT,ofile,EerrorID,Eerror
                        NCDF_VARPUT,ofile,PerrorID,Perror
                        NCDF_VARPUT,ofile,ERMqID,eRMQ
                        NCDF_VARPUT,ofile,EPLqID,ePLQ
                        NCDF_VARPUT,ofile,EEEqID,eEEQ
                        NCDF_VARPUT,ofile,EDMqID,eDMQ
                        NCDF_VARPUT,ofile,PRMqID,pRMQ
                        NCDF_VARPUT,ofile,PPLqID,pPLQ
                        NCDF_VARPUT,ofile,PEEqID,pEEQ
                        NCDF_VARPUT,ofile,PDMqID,pDMQ
                        NCDF_VARPUT,ofile,WID,W
                        NCDF_VARPUT,ofile,WpID,Wp
                        NCDF_VARPUT,ofile,pitchID,pitch
                        NCDF_VARPUT,ofile,BLC_FootID,bfofl
                        NCDF_VARPUT,ofile,BLC_SatID,blocal
                        NCDF_Varput,ofile,BLC_angleID, blc_angle

;stop
                        NCDF_CLOSE,ofile
end

;----------------------------------------------------------------------------------

function advanceDate,cdoy
	cdoy = cdoy + 1.
       	CALDAT,cdoy,dmo,ddy,dyr
	date = [string(dyr,FORMAT='(I4.4)'),string(dmo,FORMAT='(I2.2)'),string(ddy,FORMAT='(I2.2)')]
	return, date
end

;-----------------------------------------------------------------------

pro ammendNCfiles_BLC_00,indSat,date,enddate


  ;this code is to update the NC file output to include pitch angles, which are kind of important.
  ;angle = 0 ;look at 0 degree telescopes
  sat = 'POES' ;only satellite this program works with, so don't change it.
  angle = 0
  indSat = 'n16'
  ;set up individual satellite data
  if indSat eq 'n15' then begin
    ;date = ['1998','07','01'] ;year, month, day ;n15
    ;date = ['1998','07','02'] ;year, month, day ;n15
    ;endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
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
    date = ['2017','11','20'] ;year, month, day ;n18
    endDate = ['2017','11','21'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA18'
  endif
  if indSat eq 'n19' then begin
    date = ['2017','11','20'] ;year, month, day ;n19
    endDate = ['2017','11','21'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA19'
  endif
  if indSat eq 'm01' then begin
    date = ['2012','10','03'] ;year, month, day ;m02
    endDate = ['2018','12','31'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP01'
  endif
  if indSat eq 'm02' then begin
    date = ['2015','06','16'] ; actual date
    endDate = ['2015','12','31'] ; actual end date
    fullSatName = 'METOP02'
  endif
  date = ['2003','11','01'] ;year, month, day ;n16
  endDate = ['2004','07','01'] ;end date for program ;n15,16,18,19,m02
  ;fullSatName = 'NOAA16'
  
  fdoydate = JULDAY(fix(date[1]),fix(date[2]),fix(date[0]))
  enddoydate = JULDAY(fix(endDate[1]),fix(endDate[2]),fix(endDate[0]))


  while fdoydate le enddoydate do begin ;loop over days
    print, date
    ;check validity of date
    ;ncfile = '/export/home/pecked/MEPED/data/'+'POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
     ;ncfile = 'Z:\MEE_Data\Processed_MEPED_Data\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
     ncfile = 'Z:\MEE_Data\Processed_MEPED_Data\MPE_v3\Paper5\POES_combinedSpectrum_v4_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
     ;ncfile = 'Z:\MEE_Data\Allison_Data\POES_combinedSpectrum_v2_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
    exists = FILE_TEST(ncfile)
    if exists then begin
      ;open cdf data and get variables from each day using pcounts code.
      ;temp = pCounts(angle,'Electron',sat,date,indSat)

      ;read in txt data of BLC
      ;blcFile = '/export/home/pecked/MEPED/data/raw/full/nc/loss_'+indSat+date[0]+date[1]+date[2]+'.nc'
      blcFile = 'Z:\Pettit\PECK_Directory\winl26\raw\full\nc\loss_'+indSat+date[0]+date[1]+date[2]+'.nc'
      blc_exists = FILE_TEST(blcFile)
      ;stop
      if blc_exists then begin
        ;read in BLC angle from raw nc file
        blc_ncid = NCDF_OPEN(blcFile)
        varid = ncdf_varid(blc_ncid,'alpha_BLC')
        ncdf_varget,blc_ncid,varid,BLC
        NCDF_CLOSE,blc_ncid

        ;average every other time to get 16 second version of BLC
        BLC_angle = dblarr(N_ELEMENTS(BLC)/2.)
        BLC_angle = (BLC[0:-2:2]+BLC[1:-1:2])/2. ;average in 16 sec instead of indiv 8 seconds. This is fancy code!
        ;stop
                          ;open ncdf file and add in new pitch angle variable
                          BLCATTR1 = 'Bounce Loss Cone Pitch Angle'
                          BLCATTR2 = 'degrees 0 to 90'
                          ncid = NCDF_OPEN(ncfile,/WRITE)
                          NCDF_CONTROL,ncid,/REDEF
                          timedimID = 7 ;time is the last dimension and should be 7 unless someone has tampered with the original ncdf write code.
                          BLCID2 = NCDF_VARDEF(ncid,'BLC_angle',timedimID,/DOUBLE)
                          NCDF_ATTPUT,ncid,BLCID2,'Full_Name',BLCATTR1,/CHAR
                          NCDF_ATTPUT,ncid,BLCID2,'Units',BLCATTR2,/CHAR
                          NCDF_CONTROL,ncid,/ENDEF
                          NCDF_VARPUT,ncid,BLCID2,BLC_angle
                          NCDF_CLOSE,ncid
        stop
      endif else print, 'No Data this date (txt).'

    endif else print, 'No Data this date (nc).'
    ;advance to the next day
    date = advanceDate(fdoydate)
  endwhile
  
end

;-----------------------------------------------------------------------
 
pro ammendNCfiles_magField_00,indSat
  ;this code is to update the NC file output to include pitch angles, which are kind of important.
  ;angle = 0 ;look at 0 degree telescopes
  sat = 'POES' ;only satellite this program works with, so don't change it.
  angle = 0
  indSat = 'n19'
  
  ;set up individual satellite data
  
  if indSat eq 'n15' then begin
    date = ['2003','01','01']
    endDate = ['2003','10','01']
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
    date = ['2017','11','20'] ;year, month, day ;n18
    endDate = ['2017','11','21'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA18'
  endif
  if indSat eq 'n19' then begin
    date = ['2017','11','20'] ;year, month, day ;n19
    endDate = ['2017','11','21'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'NOAA19'
  endif
  if indSat eq 'm01' then begin
    date = [['2012','10','03']] ;year, month, day ;m02
    endDate = ['2018','12','31'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP01'
  endif
  if indSat eq 'm02' then begin
    ;date = ['2006','12','03'] ;year, month, day ;m02
    ;endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP02'
  endif
  if indSat eq 'm03' then begin
    ;date = ['2006','12','03'] ;year, month, day ;m02
    ;endDate = ['2014','05','11'] ;end date for program ;n15,16,18,19,m02
    fullSatName = 'METOP03'
  endif
  date = ['2004','04','02'] ;year, month, day ;n16
  endDate = ['2004','07','01'] ;end date for program ;n15,16,18,19,m02
  ;fullSatName = 'NOAA15'
  
  fdoydate = JULDAY(fix(date[1]),fix(date[2]),fix(date[0]))
  enddoydate = JULDAY(fix(endDate[1]),fix(endDate[2]),fix(endDate[0]))


  while fdoydate le enddoydate do begin ;loop over days
    print, date
    ;check validity of date
    ;ncfile = '/export/home/pecked/MEPED/data/'+'POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
    ;ncfile = 'Z:\MEE_Data\Processed_MEPED_Data\'+fullSatName+'\POES_combinedSpectrum_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
    ncfile = 'Z:\MEE_Data\Processed_MEPED_Data\MPE_v3\POES_combinedSpectrum_v4_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
    ;ncfile = 'Z:\MEE_Data\Allison_Data\POES_combinedSpectrum_v2_'+indSat+'_00_'+date[0]+date[1]+date[2]+'.nc'
    exists = FILE_TEST(ncfile)
    if exists then begin
      ;open cdf data and get variables from each day using pcounts code.
      ;temp = pCounts(angle,'Electron',sat,date,indSat)

      ;read in txt data of BLC
      ;blcFile = '/export/home/pecked/MEPED/data/raw/full/nc/loss_'+indSat+date[0]+date[1]+date[2]+'.nc'
      blcFile = 'Z:\Pettit\PECK_Directory\winl26\raw\full\nc\loss_'+indSat+date[0]+date[1]+date[2]+'.nc'
      blc_exists = FILE_TEST(blcFile)
      ;stop
      if blc_exists then begin
        ;read in BLC angle from raw nc file
        blc_ncid = NCDF_OPEN(blcFile)
        varid = ncdf_varid(blc_ncid,'Bfofl')
        ncdf_varget,blc_ncid,varid,Bfofl
        varid = ncdf_varid(blc_ncid,'Blocal')
        ncdf_varget,blc_ncid,varid,Blocal
        NCDF_CLOSE,blc_ncid

        ;average every other time to get 16 second version of BLC
        B0 = dblarr(N_ELEMENTS(Bfofl)/2.)
        Bs = dblarr(N_ELEMENTS(Blocal)/2.)
        B0 = (Bfofl[0:-2:2]+Bfofl[1:-1:2])/2. ;average in 16 sec instead of indiv 8 seconds. This is fancy code!
        Bs = (Blocal[0:-2:2]+Blocal[1:-1:2])/2. ;average in 16 sec instead of indiv 8 seconds. This is fancy code!

        ;open ncdf file and add in new pitch angle variable
        B0ATTR1 = 'Magnetic Field Strength at 120km following the foot of the field line'
        B0ATTR2 = 'Telsa?'
        BsATTR1 = 'Magnetic Field Strength at satellite measurement location'
        BsATTR2 = 'Telsa?'
        ncid = NCDF_OPEN(ncfile,/WRITE)
        NCDF_CONTROL,ncid,/REDEF
        timedimID = 7 ;time is the last dimension and should be 7 unless someone has tampered with the original ncdf write code.
        B0ID = NCDF_VARDEF(ncid,'Bfofl',timedimID,/DOUBLE)
        NCDF_ATTPUT,ncid,B0ID,'Full_Name',B0ATTR1,/CHAR
        NCDF_ATTPUT,ncid,B0ID,'Units',B0ATTR2,/CHAR
        BsID = NCDF_VARDEF(ncid,'Blocal',timedimID,/DOUBLE)
        NCDF_ATTPUT,ncid,BsID,'Full_Name',BsATTR1,/CHAR
        NCDF_ATTPUT,ncid,BsID,'Units',BsATTR2,/CHAR
        NCDF_CONTROL,ncid,/ENDEF
        NCDF_VARPUT,ncid,B0ID,B0
        NCDF_VARPUT,ncid,BsID,Bs
        NCDF_CLOSE,ncid
      endif else print, 'No Data this date (txt).'
    endif else print, 'No Data this date (nc).'
    ;advance to the next day
    date = advanceDate(fdoydate)
  endwhile

end

;-----------------------------------------------------------------------

function pCounts,degree,type,sat,date
  ;this is a recreation of a function written by Dominic Fuller-Rowell in MATLAB.
  ;This is the exact same function, except in IDL.
  filename = strlowcase(sat)+'_m02_'+string(date[0],format='(I4.4)')+string(date[1],format='(I2.2)')+string(date[2],format='(I2.2)')+'.cdf'

  if type eq 'Proton' then begin
    if degree eq 0 then begin
      elem = 11
      elem1 = 11
      elem2 = 5
    endif
    if degree eq 90 then begin
      elem1 = 11
      elem = 12
      elem2 = 6
    endif
    if degree ne 0 and degree ne 90 then begin
      print, 'Please choose either 0 or 90 degree telescope'
      return,0
    endif
  endif

  if type eq 'Electron' then begin
    if degree eq 0 then begin
      elem = 13
      elem1 = 11
      elem2 = 5
    endif
    if degree eq 90 then begin
      elem = 14
      elem1 = 12
      elem2 = 6
    endif
  endif

  if type eq 'Protron_not' then begin

    if degree eq 0 then begin
      elem = 11
      elem2 = 5
    endif

    if degree eq 90 then begin
      elem = 12
      elem2 = 6
    endif

    if degree ne 0 and degree ne 90 then begin
      print, 'Please choose either 0 or 90 degree telescope'
      return,0
    endif
  endif

  if type eq 'Electron_not' then begin
    if degree eq 0 then begin
      elem = 13
      elem2 = 5
    endif
    if degree eq 90 then begin
      elem = 14
      elem2 = 6
    endif
    if degree ne 0 and degree ne 90 then begin
      print, 'Please choose either 0 or 90 degree telescope'
      return,0
    endif
  endif
  ;if type ne 'Proton' and type ne 'Electron' then begin
  ;  print, 'Please choose either Electron or Proton response'
  ;  return,0
  ;endif
  cdf_exists = FILE_TEST('D:\POES_Data\'+filename)
  if cdf_exists then begin
  cdfid = CDF_OPEN('D:\POES_Data\'+filename,/READONLY)
  CDF_CONTROL,cdfid,variable=elem,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,elem,particleCount,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=elem1,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,elem1,pparticleCount,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=1,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,1,geogLL,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=0,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,0,time,/ZVARIABLE,REC_COUNT=(info.MAXREC+1) ;note that you need to use cdf_epoch,time[i],yr,mn,dy,hr,min,sec,milli,/BREAKDOWN_EPOCH to make sense of this variable.
  CDF_CONTROL,cdfid,variable=3,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,3,lValue,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=2,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,2,foflLL,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=4,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,4,MLT,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  CDF_CONTROL,cdfid,variable=elem2,/zvariable,get_var_info=info
  CDF_VARGET,cdfid,elem2,pitch,/ZVARIABLE,REC_COUNT=(info.MAXREC+1)
  MLT = reform(MLT)
  lvalue = reform(lvalue)
  time = reform(time)
  pitch = reform(pitch)
  
 ; stop
  save, pparticleCount, particleCount, geogLL, time, MLT, foflLL, pitch, lValue, filename='D:\POES_Data\MPE_Software\0_Degree_Detector\Pcounts_temp.sav'
  ;stop
  return,CREATE_STRUCT('particleCount',REFORM(particleCount),'geogLL',REFORM(geogLL),'time',REFORM(time),'lValue',REFORM(lvalue),'foflLL',REFORM(foflLL),'MLT',REFORM(MLT),'pitch',REFORM(pitch))
  endif else particleCount = -999
  return, CREATE_STRUCT('particleCount', particleCount) 
end


function pCounts2,degree,type,sat,date

  degree='0'
  sat = 'POES'
  type = 'Electron'
  satname = 'METOP03'
  indsat = 'm03'

  ifilename = 'D:\POES_Data\New_POES_Data\'+satname+'\'+strlowcase(sat)+'_'+indsat+'_'+string(date[0],format='(I4.4)')+string(date[1],format='(I2.2)')+string(date[2],format='(I2.2)')+'_raw.nc'
  ;filename = 'Z:\MEE_Data\POES_Data\New_POES_Data\poes_n17_20120201_raw.nc'
  
  cdf_exists = FILE_TEST(ifilename)
  if cdf_exists then begin
  cdfid = ncdf_open(ifilename)
  ;  stop
  ; ---------------------------------------------------------------------
  ; Read data
  ; ---------------------------------------------------------------------
    varid = ncdf_varid(cdfid, 'mep_ele_tel0_cps_e1') ;electron telescope channel 1
    ncdf_varget, cdfid, varid, mep_ele_tel0_cps_e1
    varid = ncdf_varid(cdfid, 'mep_ele_tel0_cps_e2') ;electron telescope channel 2
    ncdf_varget, cdfid, varid, mep_ele_tel0_cps_e2
    varid = ncdf_varid(cdfid, 'mep_ele_tel0_cps_e3') ;electron telescope channel 3
    ncdf_varget, cdfid, varid, mep_ele_tel0_cps_e3
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p1') ;proton telescope channel 1
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p1
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p2') ;proton telescope channel 2
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p2
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p3') ;proton telescope channel 3
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p3
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p4') ;proton telescope channel 4
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p4
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p5') ;proton telescope channel 5
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p5
    varid = ncdf_varid(cdfid, 'mep_pro_tel0_cps_p6') ;proton telescope channel 6
    ncdf_varget, cdfid, varid, mep_pro_tel0_cps_p6
  varid = ncdf_varid(cdfid, 'time')                ; time in milliseconds since 1970
  ncdf_varget, cdfid, varid, time
  varid = ncdf_varid(cdfid, 'day')                 ; day of year  DDD
  ncdf_varget, cdfid, varid, day
  varid = ncdf_varid(cdfid, 'year')                ; year   YYYY
  ncdf_varget, cdfid, varid, year
  varid = ncdf_varid(cdfid, 'lat')  ; Latitude of satellite = geogLat
  ncdf_varget, cdfid, varid, lat
  varid = ncdf_varid(cdfid, 'lon')  ; Longitude of satellite = geogLon
  ncdf_varget, cdfid, varid, lon
  varid = ncdf_varid(cdfid, 'time')                         ; time
  ncdf_varget, cdfid, varid, time
  ncdf_close,cdfid

  newtime = time-time[0]
  newtime = newtime/3600000.

  ifilename2 = 'D:\POES_Data\New_POES_Data\'+satname+'\'+strlowcase(sat)+'_'+indsat+'_'+string(date[0],format='(I4.4)')+string(date[1],format='(I2.2)')+string(date[2],format='(I2.2)')+'_proc.nc'
  ;filename2 = 'Z:\MEE_Data\POES_Data\New_POES_Data\poes_n17_20120201_proc.nc'

  cdfid2 = ncdf_open(ifilename2)
  varid = ncdf_varid(cdfid2, 'geod_lat_foot')                ; geodetic latitude (foot of field line)
  ncdf_varget, cdfid2, varid, geod_lat_foot
  varid = ncdf_varid(cdfid2, 'geod_lon_foot')                ; geodetic longitude (foot of field line)
  ncdf_varget, cdfid2, varid, geod_lon_foot
  varid = ncdf_varid(cdfid2, 'MLT')     ; MLT
  ncdf_varget, cdfid2, varid, MLT
  varid = ncdf_varid(cdfid2, 'L_IGRF')  ; L value
  ncdf_varget, cdfid2, varid, L_IGRF
  varid = ncdf_varid(cdfid2, 'meped_alpha_0_sat')                ; pitch angle of the satellite
  ncdf_varget, cdfid2, varid, meped_alpha_0_sat
  varid = ncdf_varid(cdfid2, 'mep_IFC_on')                ; pitch angle of the satellite
  ncdf_varget, cdfid2, varid, mep_IFC_on
  varid = ncdf_varid(cdfid2, 'lat')                          ; Satellite Latitude
  ncdf_varget, cdfid2, varid, lat
  varid = ncdf_varid(cdfid2, 'lon')                          ; Satellite Longitude
  ncdf_varget, cdfid2, varid, lon
  varid = ncdf_varid(cdfid2, 'Btot_sat')                          ; magnetic field at satellite
  ncdf_varget, cdfid2, varid, Btot_sat
  varid = ncdf_varid(cdfid2, 'Btot_foot')                          ; magnetic field at footprint of satellite
  ncdf_varget, cdfid2, varid, Btot_foot
  ncdf_close,cdfid

  if N_ELEMENTS(lat) eq N_ELEMENTS(time) then begin

  ; average the data into 16secs
  ;stop
  end_time = N_ELEMENTS(time)/8.
  end_time = round(end_time)

  E1_Channel_avg = fltarr(end_time)
  E2_Channel_avg = fltarr(end_time)
  E3_Channel_avg = fltarr(end_time)
  P1_Channel_avg = fltarr(end_time)
  P2_Channel_avg = fltarr(end_time)
  P3_Channel_avg = fltarr(end_time)
  P4_Channel_avg = fltarr(end_time)
  P5_Channel_avg = fltarr(end_time)
  P6_Channel_avg = fltarr(end_time)
  l_value = fltarr(end_time)
  l_value_avg = fltarr(end_time)
  pitch_angle = fltarr(end_time)
  time_avg = fltarr(end_time)
  index_array = []
  MLT_avg = fltarr(end_time)
  lat_avg = fltarr(End_time)
  lon_avg = fltarr(end_time)
  folat_avg = fltarr(end_time)
  folon_avg = fltarr(end_time)
  msec_time = fltarr(end_time)
  B_Total_sat_avg = fltarr(end_time)
  B_Total_foot_avg = fltarr(end_time)
  BLC_Alpha = fltarr(end_time)

  for i = 0, end_time-1 do begin

    if i eq 0 then index1 = 0 ;3, 11
    if i eq 0 then index2 = 7 ;10, 18
    if i eq 0 then index3 = 4 ;15

 ;   if i eq 1 then index1 = 19 ;11
 ;   if i eq 1 then index2 = 26 ;18
 ;   if i eq 1 then index3 = 23

    if i ge 1 then index1 = index1 + 8.
    if i ge 1 then index2 = index2 + 8.
    if i ge 1 then index3 = index3 + 8.

    if index1 ge N_ELEMENTS(time) then index1 = N_ELEMENTS(time)-1
    if index2 ge N_ELEMENTS(time) then index2 = N_ELEMENTS(time)-1
    if index3 ge N_ELEMENTS(time) then index3 = N_ELEMENTS(time)-1

    Channel_1_Temp = mep_ele_tel0_cps_e1[index1:index2]
    Channel_2_Temp = mep_ele_tel0_cps_e2[index1:index2]
    Channel_3_Temp = mep_ele_tel0_cps_e3[index1:index2]
    PChannel_1_Temp = mep_pro_tel0_cps_p1[index1:index2]
    PChannel_2_Temp = mep_pro_tel0_cps_p2[index1:index2]
    PChannel_3_Temp = mep_pro_tel0_cps_p3[index1:index2]
    PChannel_4_Temp = mep_pro_tel0_cps_p4[index1:index2]
    PChannel_5_Temp = mep_pro_tel0_cps_p5[index1:index2]
    PChannel_6_Temp = mep_pro_tel0_cps_p6[index1:index2]
;    Error_Flag_temp = mep_ifc_on[index1:index2]
    B_Total_sat_temp = Btot_sat[index1:index2]
    B_Total_foot_temp = Btot_foot[index1:index2]
 ;   stop
 ;   error_index = where(error_flag_temp eq 0)
 ;   channel_1_temp2 = channel_1_temp
 ;   channel_1_temp2[error_index] = 0./0.
;stop
    E1_Channel_avg[i] =  mean(Channel_1_Temp, /naN, /double)
    E2_Channel_avg[i] =  mean(Channel_2_Temp, /naN, /double)
    E3_Channel_avg[i] =  mean(Channel_3_Temp, /naN, /double)
    P1_Channel_avg[i] =  mean(PChannel_1_Temp, /naN, /double)
    P2_Channel_avg[i] =  mean(PChannel_2_Temp, /naN, /double)
    P3_Channel_avg[i] =  mean(PChannel_3_Temp, /naN, /double)
    P4_Channel_avg[i] =  mean(PChannel_4_Temp, /naN, /double)
    P5_Channel_avg[i] =  mean(PChannel_5_Temp, /naN, /double)
    P6_Channel_avg[i] =  mean(PChannel_6_Temp, /naN, /double)
    B_Total_sat_avg[i]=  mean(B_Total_sat_temp, /NaN, /double)
    B_Total_foot_avg[i]=  mean(B_Total_foot_temp, /NaN, /double)

    lat_avg[i]       =  lat[index3]
    lon_avg[i]       =  lon[index3]
    folat_avg[i]       =  geod_lat_foot[index3]
    folon_avg[i]       =  geod_lon_foot[index3]

    pitch_angle[i]    =  meped_alpha_0_sat[index3]
    l_value_avg[i]    =  mean(L_IGRF[index1:index2], /Nan, /double)
    l_value[i]        =  L_IGRF[index3]
    mlt_avg[i]        =  MLT[index3]
    time_avg[i]       =  newtime[index3]
    msec_time[i]      =  time[index3]

    if l_value[i] ge 20. then l_value[i] = -999.

    index_temp = [index1, index2]
    index_array = [[index_array], [index_temp]]

  endfor

  ; Process the data into expected arrays

  particleCount = [[e1_channel_avg], [e2_channel_avg], [e3_channel_avg]]
  particleCount = transpose(particleCount)
  pparticleCount = [[p1_channel_avg], [p2_channel_avg], [p3_channel_avg], [p4_channel_avg], [p5_channel_avg], [p6_channel_avg]]
  pparticleCount = transpose(pparticleCount)

  geogLL = [[lat_avg], [lon_avg]]
  geogLL = transpose(geogLL)
  foflLL = [[folat_avg], [folon_avg]]
  foflLL = transpose(foflLL)
  pitch = pitch_angle
  lValue = l_value
  MLT = mlt_avg
  time = time_avg
  Bfofl = B_TOTAL_FOOT_AVG
  Blocal = B_Total_Sat_AVG

  BLC_Alpha_temp = sqrt((B_TOTAL_SAT_AVG/B_TOTAL_FOOT_AVG))
  BLC_alpha = asin(BLC_ALPHA_Temp)

  BLC_Angle = blc_alpha*180./!pi
  ;time[0] = time[0]/2.
  
  ; stop
  save, BLC_Angle, Bfofl, Blocal, pparticleCount, particleCount, msec_time, geogLL, time, MLT, foflLL, pitch, lValue, filename='D:\POES_Data\MPE_Software\0_Degree_Detector\Pcounts_temp.sav'
  return,CREATE_STRUCT('Blc_angle',reform(blc_angle),'Bfofl',REFORM(Bfofl),'Blocal',REFORM(Blocal),'particleCount',REFORM(particleCount),'geogLL',REFORM(geogLL),'time',REFORM(time),'lValue',REFORM(lvalue),'foflLL',REFORM(foflLL),'MLT',REFORM(MLT),'pitch',REFORM(pitch))
  endif else particleCount = -999
  endif else particleCount = -999
  return,CREATE_STRUCT('particleCount',REFORM(particleCount))
end
