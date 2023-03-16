ifiles_range = file_search('Z:\MEE_Data\Processed_MEPED_Data\POES_combinedSpectrum_m01_00_20151003.nc')


cdfid = ncdf_open(ifiles_range)
;varid = ncdf_varid(cdfid2, 'BLC_angle')
;ncdf_varget, cdfid2, varid, BLC_angle
varid = ncdf_varid(cdfid, 'time')
ncdf_varget, cdfid, varid, time
varid = ncdf_varid(cdfid, 'Ecounts')
ncdf_varget, cdfid, varid, Ecounts
varid = ncdf_varid(cdfid, 'EDMq')
ncdf_varget, cdfid, varid, EDMq
varid = ncdf_varid(cdfid, 'EEEq')
ncdf_varget, cdfid, varid, EEEQ
varid = ncdf_varid(cdfid, 'Eerror')
ncdf_varget, cdfid, varid, Error
varid = ncdf_varid(cdfid, 'energy')
ncdf_varget, cdfid, varid, energy
varid = ncdf_varid(cdfid, 'EOcounts')
ncdf_varget, cdfid, varid, EOcounts
varid = ncdf_varid(cdfid, 'EOcounts_corrected')
ncdf_varget, cdfid, varid, EOcounts_corrected
varid = ncdf_varid(cdfid, 'EOcounts_Lambda')
ncdf_varget, cdfid, varid, EOcounts_Lambda
varid = ncdf_varid(cdfid, 'EPLq')
ncdf_varget, cdfid, varid, EPLq
varid = ncdf_varid(cdfid, 'foflLat')
ncdf_varget, cdfid, varid, foflLat
varid = ncdf_varid(cdfid, 'foflLon')
ncdf_varget, cdfid, varid, foflLon
varid = ncdf_varid(cdfid, 'geogLat')
ncdf_varget, cdfid, varid, geogLat
varid = ncdf_varid(cdfid, 'geogLon')
ncdf_varget, cdfid, varid, geogLon
varid = ncdf_varid(cdfid, 'lValue')
ncdf_varget, cdfid, varid, lValue
varid = ncdf_varid(cdfid, 'MLT')
ncdf_varget, cdfid, varid, MLT
varid = ncdf_varid(cdfid, 'Pcounts')
ncdf_varget, cdfid, varid, Pcounts
varid = ncdf_varid(cdfid, 'PDMq')
ncdf_varget, cdfid, varid, PDMq
varid = ncdf_varid(cdfid, 'PEEq')
ncdf_varget, cdfid, varid, PEEq
varid = ncdf_varid(cdfid, 'Perror')
ncdf_varget, cdfid, varid, Perror
varid = ncdf_varid(cdfid, 'pitch')
ncdf_varget, cdfid, varid, pitch
varid = ncdf_varid(cdfid, 'POcounts')
ncdf_varget, cdfid, varid, POcounts
varid = ncdf_varid(cdfid, 'POcounts_Lambda')
ncdf_varget, cdfid, varid, POcounts_Lambda
varid = ncdf_varid(cdfid, 'PPLq')
ncdf_varget, cdfid, varid, PPLq
varid = ncdf_varid(cdfid, 'PRMq')
ncdf_varget, cdfid, varid, PRMq
varid = ncdf_varid(cdfid, 'time')
ncdf_varget, cdfid, varid, time
varid = ncdf_varid(cdfid, 'W')
ncdf_varget, cdfid, varid, W
varid = ncdf_varid(cdfid, 'Wp')
ncdf_varget, cdfid, varid, Wp


stop

ifiles_range2 = file_search('Z:\MEE_Data\Processed_MEPED_Data\POES_combinedSpectrum_m01_00_v1_20151003.nc')


cdfid2 = ncdf_open(ifiles_range2)
;varid = ncdf_varid(cdfid2, 'BLC_angle')
;ncdf_varget, cdfid2, varid, BLC_angle2
;varid = ncdf_varid(cdfid2, 'Bfofl')
;ncdf_varget, cdfid2, varid, Bfofl
;varid = ncdf_varid(cdfid2, 'Blocal')
;ncdf_varget, cdfid2, varid, Blocal
varid = ncdf_varid(cdfid, 'time')
ncdf_varget, cdfid, varid, time
varid = ncdf_varid(cdfid2, 'Ecounts')
ncdf_varget, cdfid2, varid, Ecounts2
varid = ncdf_varid(cdfid2, 'EDMq')
ncdf_varget, cdfid2, varid, EDMq2
varid = ncdf_varid(cdfid2, 'EEEq')
ncdf_varget, cdfid2, varid, EEEQ2
varid = ncdf_varid(cdfid2, 'Eerror')
ncdf_varget, cdfid2, varid, Error2
varid = ncdf_varid(cdfid2, 'energy')
ncdf_varget, cdfid2, varid, energy2
varid = ncdf_varid(cdfid2, 'EOcounts')
ncdf_varget, cdfid2, varid, EOcounts2
varid = ncdf_varid(cdfid2, 'EOcounts_corrected')
ncdf_varget, cdfid2, varid, EOcounts_corrected2
varid = ncdf_varid(cdfid2, 'EOcounts_Lambda')
ncdf_varget, cdfid2, varid, EOcounts_Lambda2
varid = ncdf_varid(cdfid2, 'EPLq')
ncdf_varget, cdfid2, varid, EPLq2
varid = ncdf_varid(cdfid2, 'foflLat')
ncdf_varget, cdfid2, varid, foflLat2
varid = ncdf_varid(cdfid2, 'foflLon')
ncdf_varget, cdfid2, varid, foflLon2
varid = ncdf_varid(cdfid2, 'geogLat')
ncdf_varget, cdfid2, varid, geogLat2
varid = ncdf_varid(cdfid2, 'geogLon')
ncdf_varget, cdfid2, varid, geogLon2
varid = ncdf_varid(cdfid2, 'lValue')
ncdf_varget, cdfid2, varid, lValue2
varid = ncdf_varid(cdfid2, 'MLT')
ncdf_varget, cdfid2, varid, MLT2
varid = ncdf_varid(cdfid2, 'Pcounts')
ncdf_varget, cdfid2, varid, Pcounts2
varid = ncdf_varid(cdfid2, 'PDMq')
ncdf_varget, cdfid2, varid, PDMq2
varid = ncdf_varid(cdfid2, 'PEEq')
ncdf_varget, cdfid2, varid, PEEq2
varid = ncdf_varid(cdfid2, 'Perror')
ncdf_varget, cdfid2, varid, Perror2
varid = ncdf_varid(cdfid2, 'pitch')
ncdf_varget, cdfid2, varid, pitch2
varid = ncdf_varid(cdfid2, 'POcounts')
ncdf_varget, cdfid2, varid, POcounts2
varid = ncdf_varid(cdfid2, 'POcounts_Lambda')
ncdf_varget, cdfid2, varid, POcounts_Lambda2
varid = ncdf_varid(cdfid2, 'PPLq')
ncdf_varget, cdfid2, varid, PPLq2
varid = ncdf_varid(cdfid2, 'PRMq')
ncdf_varget, cdfid2, varid, PRMq2
varid = ncdf_varid(cdfid2, 'time')
ncdf_varget, cdfid2, varid, time2
varid = ncdf_varid(cdfid2, 'W')
ncdf_varget, cdfid2, varid, W2
varid = ncdf_varid(cdfid2, 'Wp')
ncdf_varget, cdfid2, varid, Wp2

stop

for i = 0, 5399 do begin
  
  if eocounts[0,i] gt 100 then begin
    if eocounts2[0, i] gt 100 then begin
      
      if eocounts[2,i] ne eocounts2[2,i] then begin
  
        ;a = plot(energy, eocounts[*, i], color='blue', /ylog)
        ;b = plot(energy, eocounts2[*, i], color='red', /overplot, linestyle='dash')
        print, eocounts[*, i]
        print, eocounts_corrected[*,i]
        print, eocounts_corrected2[*,i]    
        print, pcounts[*, i]
      
        a = plot(energy, ecounts[*,i], color='blue', /ylog, /xlog)
        b = plot(energy, ecounts2[*,i], color='red', /overplot, linestyle='dash')
        stop
      endif 
    endif
  endif
 
endfor

stop

ifiles_range3 = file_search('Z:\MEE_Data\POES_Data\New_POES_Data\poes_n15_20131120_proc.nc')
cdfid3 = ncdf_open(ifiles_range3)

varid = ncdf_varid(cdfid3, 'Btot_foot')
ncdf_varget, cdfid3, varid, Btot_foot
varid = ncdf_varid(cdfid3, 'Btot_sat')
ncdf_varget, cdfid3, varid, Btot_sat


stop
end
