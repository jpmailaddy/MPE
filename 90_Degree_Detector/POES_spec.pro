;POES_spec.pro
;created by, Ethan D. Peck
;This code contains some POES specific codes (hence POES_spec). 

function validSats, date
;this function will check the date against available satellites for that date. That way the program will search for it.
	sats = ['n15','n16','n17','n18','n19','m02'] ;names of all POES satellites
	sdate = [''] ;put start dates for each satellite here
	edate = [''] ;put end dates for each satellite here ;will put Jan 1st, 2200 for satellites that are still providing data. I assume they will be dead by then... update this as they stop working.
	return, 0 ;for now I am not using this function. Eventually I want this to return an array of satellite names from sats variable that are applicable to the given date.

end


function pCounts,degree,type,sat,date
  ;this is a recreation of a function written by Dominic Fuller-Rowell in MATLAB.
  ;This is the exact same function, except in IDL.
  filename = strlowcase(sat)+'_n15_'+string(date[0],format='(I4.4)')+string(date[1],format='(I2.2)')+string(date[2],format='(I2.2)')+'.cdf'

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
 
  if type eq 'Electron' then begin
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

  cdfid = CDF_OPEN('Z:\Pettit\PECK_Directory\winl26\raw\'+filename,/READONLY)
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
  ;stop
  save, pparticleCount, particleCount, geogLL, time, MLT, foflLL, pitch, lValue, filename='Z:\Pettit\RAISE_MEE_Simulations\PECK_MAPS\Final\Pcounts_temp.sav'
  ;stop
  return,CREATE_STRUCT('particleCount',REFORM(particleCount),'geogLL',REFORM(geogLL),'time',REFORM(time),'lValue',REFORM(lvalue),'foflLL',REFORM(foflLL),'MLT',REFORM(MLT),'pitch',REFORM(pitch))
end

function readResponse, TelescopeType,ResponseType,minEnergy,maxProtonEnergy,maxElectronEnergy,interp
;this is a recreation of a function written by Dominic Fuller-Rowell in MATLAB.
;This is the exact same function, except in IDL.
;Function to read and log linearly interpolate the response functions from CSV files in ./data
;I believe this is just the G factor Tables from Yando et al. 2011.
        ;if TelescopeType eq 'Proton' and ResponseType eq 'Proton' then debugmode = 1 else debugmode = 0
        sat = 'POES'
        ;I will differ from original matlab code here to take advantage of CSV reading function built into idl
        ;data = READ_CSV('./data/'+sat+'_'+TelescopeType+'Telescope_'+ResponseType+'Response.csv')
        data = READ_CSV('D:\POES_Data\MPE_Software\Resources\'+sat+'_'+TelescopeType+'Telescope_'+ResponseType+'Response.csv')
        energy = data.field1
        response = dblarr(N_ELEMENTS(energy),N_TAGS(data)-2)
        for i = 1,N_TAGS(data)-2 do response[*,i-1]=data.(i) ;loop through tags of data structure (skipping first tag as that is the energy field.
        ;Calculate logarithmic energy axis from min and max energies chosen in the run.m
        edges = logspace(alog10(minEnergy),alog10(maxProtonEnergy),interp+1)
        midpoints = dblarr(interp)
        deltaE = dblarr(interp)
        response[where(response eq 0,/null)] = 1.e-18 ;Do not allow 0 values in response function ;GFactor or GlnbFactor
        ;response[where(response eq 0 OR FINITE(response) eq 0,/null)] = 1.e-18 ;Do not allow 0 values in response function ;GlnFactor
        ;Calculate midpoints and deltaE for adjacent energies.
        for i = 0, interp-1 do begin
                midpoints[i] = sqrt(edges[i]*edges[i+1])
                deltaE[i] = edges[i+1]-edges[i]
        endfor
        ;Choose subset of energy axis for Electron response
        if ResponseType eq 'Electron' then begin
                element = where(abs(midPoints-maxElectronEnergy) eq min(abs(midpoints-maxElectronEnergy)),/null)
                maxElectronEnergytemp = (midPoints[element])[0]
                midPoints = midPoints[0:where(long(midPoints) eq long(maxElectronEnergytemp),/null)] ;hopefully this will work in the same way as int16 in the matlab code.
        end
        ;create virtual indices for interpol method.
        JY = dblarr(N_ELEMENTS(midpoints))
        for ind = 0,N_ELEMENTS(midpoints)-1 do begin
                tind1 = double((where(energy ge midpoints[ind],/null))[0]) ;first index of tind1 is just above midpoints
                tind2 = double((where(energy le midpoints[ind],/null))[-1]) ;last index of tind2 is just below midpoints
                JY[ind] = double((midpoints[ind]-energy[tind2])/(energy[tind1]-energy[tind2]))+tind2 ;will normalize midpoints between surrounding energies and add normalization to current index to create virtual indeces. Not sure this will work the way I think it will.
        endfor
        response = BILINEAR(alog10(response),JY,indgen(N_ELEMENTS(response[0,*]))) ;GFactor
        response = 10.^response ;GFactor
        s = size(response)
        x = s[1]
        y = s[2]
        newDeltaE = dblarr(x,y)
        for i = 0,y-1 do begin ;note that in matlab code this was a for loop of y not x. So my matrix is reversed.
                newDeltaE[*,i] = deltaE[0:x-1]
        endfor
        ;if debugmode then stop 
        ;Undo setting values to non-zero
        for i = 0,x-1 do begin
                for j = 0,y-1 do begin
                        if response[i,j] lt 1e-17 then response[i,j] = 0 ;shouldnt it be less than 1e-3? ;come here again.
                endfor
        endfor
        energy = midpoints
        ;if debugmode then stop
        response = response*newDeltaE ;multiply responses by their corresponding deltaE
        ;if debugmode then stop
        response = response/100. ;convert values in cm^2*sr instead of 100*cm^2*sr.
        ;if debugmode then stop
        return,CREATE_STRUCT('response',response,'energy',energy,'deltaE',deltaE)
end

function logspace, A, B, N
;idl version of the logspace function from matlab.
;creates n value locations equally separated between positions A and B in log space.
        L = DINDGEN(N) / (N - 1.0D) * (B - A) + A
        L = 10^L
        return, L
end

function logNorm,covarElem,energyAxis,sigma
;this is a recreation of a function written by Dominic Fuller-Rowell in MATLAB.
;This is the exact same function, except in IDL.
        ;print, where(covarElem lt 0)
        numElems = N_ELEMENTS(covarElem)
        covarMat = dblarr(numElems,numElems)
        for i = 0,numElems-1 do begin
                for j = 0,numElems-1 do begin
                        exponent = (((alog(energyAxis[i])-alog(energyAxis[j]))^2.)/(-2.*sigma^2.))
                        covarMat[i,j] = (sqrt(covarElem[i]*covarElem[j])*exp(exponent))
                        if finite(covarMat[i,j]) eq 0 then stop
                endfor
        endfor
        return, covarMat
end

function simpleForwardModel, K, f, dt,E
; lambda = Hf ;eq. 2
; H = dt*G*dE = K*dE ;eq. 3
; dt = integration time (16 sec. for POES 16-sec data)
; dE is trapazoid rule energy range.
; K = H = dt*G*dE, weighting function (geometric factor * time), cm2*sr*sec
; f is the combined solution (counts/cm2/sr/keV)
; lambda is expected counts (counts)

  ;Ksize = size(K)
  ;H = dblarr(Ksize[1],Ksize[2])
  ;dE = dblarr(N_ELEMENTS(E)) ;eq. 5
  ;dE[0] = (E[1]-E[0])/2.
  ;dE[-1] = (E[-1]-E[-2])/2.
  ;dE[1:-2] = (E[2:-1]-E[0:-3])/2.
  ;for i = 0,Ksize[2]-1 do begin ;number of channels loop
  ;  H[*,i] = (K[*,i]*dE) ;eq. 3
  ;endfor 
  ;lambda = REFORM(H##f)
  lambda = REFORM(K##f)
  return, lambda

end
 

