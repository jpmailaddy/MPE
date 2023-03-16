function gf_me1_recommended, energy

; From GOESN-ENG-028, Table 3-54 (with one point from Table 3-49, averaged 0E1 and 90E1 GF at 20 keV)
; energy = electron energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 8
   
   xenergy = [20.0, 25.0, 35.0, 45.0, 60.0, 200.0, 500.0, 2000.0]*1.0d0

   gf = [8.2e-8, 3.0e-5, 1.0d-2, 1.0d-2, 1.0d-3, 1.0d-3, 2.0d-4, 2.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), alog10(xenergy), alog10(energy))
   endelse


end

;------------------------------------------------------------------------------------------

function gf_me2_recommended, energy

; From GOESN-ENG-028, Table 3-54 (with one point from Table 3.43)
; energy = electron energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 8
   
   xenergy = [40.0, 45.0, 55.0, 95.0, 110.0, 250.0, 500.0, 2000.0]*1.0d0

   gf = [1.03e-6, 3.0e-5, 1.0d-2, 1.0d-2, 1.0d-3, 1.0d-3, 4.0d-4, 4.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_me3_recommended, energy

; From GOESN-ENG-028, Table 3-54 (with four points from Table 3.43)
; energy = electron energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 9
   
   xenergy = [55.0, 60.0, 75.0, 90.0, 95.0, 105.0, 200.0, 500.0, 2000.0]*1.0d0

   gf = [5.75d-7, 1.42d-6, 2.24d-6, 1.39d-6, 1.0d-4, 1.0d-2, 1.0d-2, 1.0d-3, 1.0d-3]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_me4_recommended, energy

; From GOESN-ENG-028, Table 3-54 (with one point from Table 3.45)
; energy = electron energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 7
   
   xenergy = [114.0, 190.0, 220.0, 350.0, 500.0, 800.0, 2000.0]*1.0d0

   gf = [4.87e-4, 1.0d-3, 1.0d-2, 1.0d-2, 3.0d-3, 6.0d-3, 1.0d-2]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
;   endif else if energy ge xenergy[2] and energy le xenergy[3] then begin
;      return, 1.0d-2
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_me5_recommended, energy

; From GOESN-ENG-028, Table 3-54 (with three points from Table 3.45)
; energy = electron energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 8
   
   xenergy = [114.0, 168.0, 261.0, 350.0, 400.0, 600.0, 1000.0, 2000.0]*1.0d0

   gf = [1.0d-7, 6.94d-7, 2.99d-7, 1.0d-3, 1.0d-2, 1.0d-2, 4.0d-3, 3.0d-3]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

pro maged_gf, energy

   energy = findgen(401)*5.0d0 ; 0-2000 keV
   gf = dblarr(401, 5)
;   energy = findgen(1001)*2.0d0 ; 0-2000 keV
;   gf = dblarr(1001, 5)
   
   for i = 0, 400 do begin
      gf[i,0] = gf_me1_recommended(energy[i])
      gf[i,1] = gf_me2_recommended(energy[i])
      gf[i,2] = gf_me3_recommended(energy[i])
      gf[i,3] = gf_me4_recommended(energy[i])
      gf[i,4] = gf_me5_recommended(energy[i])
   endfor
   
   set_plot, 'win'

   colorcat12, color_scale
   !p.multi = [0,1,1]

   window, 10
   
   xrng = [20, 2000]
   yrng = [1e-7, 1e-1]
   
   plot, energy, gf[*,0], /nodata, /xlog, /ylog, xr = xrng, yr = yrng, xstyle = 1, ystyle = 1, $
      xtitle = 'Electron Energy (keV)', ytitle = 'Geometrical Factor (cm!U2!N sr)', $
      title = 'MAGED Geometrical Factors'

   for j = 0, 4 do oplot, energy, gf[*,j], color = color_scale[(2*j+1)], thick = 2

   Img = TVRD(0, TRUE = 1)
   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\bowtie\maged_gf_loglog.png'
   WRITE_PNG, PlotnamePNG, Img
   
   window, 11
   
   xrng = [10, 2000]
   yrng = [0, 1.2e-2]
   
   plot, energy, gf[*,0], /nodata, /xlog, xr = xrng, yr = yrng, xstyle = 1, ystyle = 1, $
      xtitle = 'Electron Energy (keV)', ytitle = 'Geometrical Factor (cm!U2!N sr)', $
      title = 'MAGED Geometrical Factors'

   for j = 0, 4 do oplot, energy, gf[*,j], color = color_scale[(2*j+1)], thick = 2

   Img = TVRD(0, TRUE = 1)
   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\bowtie\maged_gf_linlog.png'
   WRITE_PNG, PlotnamePNG, Img
   
;   Img = TVRD(0, TRUE = 1)
;   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\PAD inversion\weighting_function.png'
;   PlotnamePNG = 'C:\Data\GOES-R\PAD inversion\weighting_function.png'

;   WRITE_PNG, PlotnamePNG, Img

; stop
      
end

