;------------------------------------------------------------------------------------------

function gf_mp1_recommended, energy

; From GOESN-ENG-029, Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 4
   
   xenergy = [75., 85., 105., 115.]*1.0d0

   gf = [1.0d-4, 1.0d-2, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp1_measured, energy

; From GOESN-ENG-029, From GOESN-ENG-029, Table 3-29
; energy = proton energy in keV
; geometrical factor in cm^2 sr

; might be ~2x low due to beam misalignment

   energy = double(energy)
   
   nenergy = 15
   
   xenergy = [69.1, 74.2, 79.3, 85.0, 97.8, 106.9, 111.3, 116.7, $
      122.5, 157.0, 200., 218., 241., 274., 521.]*1.0d0

   gf = [7.26d-7, 8.06d-5, 1.12d-3, 4.34d-3, 5.92d-3, 4.71d-3, 4.37d-3, 6.74d-4, $
      6.21d-5, 1.48d-5, 6.48d-6, 6.41d-6, 8.01d-6, 5.72d-6, 1.31d-5]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), alog10(xenergy), alog10(energy))
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp2_recommended, energy

; From GOESN-ENG-029, Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 4
   
   xenergy = [105., 115., 165., 175.]*1.0d0

   gf = [1.0d-4, 1.0d-2, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp2_composite, energy

; From GOESN-ENG-029, Table 3-29 and Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 15
   
   xenergy = [74.2, 79.3, 85.0, 97.8, 106.9, 111.3, 116.7, 122.5, 157.0, $
      165., 175., 200., 218., 241., 274., 521.]*1.0d0

   gf = [5.19d-7, 4.70d-7, 2.35d-7, 1.46d-6, 1.7d-4, 2.06d-3, 5.59d-3, 6.13d-3, 1.62d-2, $
      1.0d-2, 1.0d-4, 2.0d-5, 1.96d-5, 1.18d-5, 8.05d-6, 3.77d-5]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp3_recommended, energy

; From GOESN-ENG-029, Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 4
   
   xenergy = [165., 175., 245., 255.]*1.0d0

   gf = [1.0d-4, 1.0d-2, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp3_composite, energy

; From GOESN-ENG-029, Table 3-29 and Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 11
   
   xenergy = [122.5, 157.0, 165., 175., $
      200., 218., 241., 245., 255., 274., 521.]*1.0d0

   gf = [1.53d-7, 1.37d-5, 1.0d-4, 1.0d-2, $
      1.84d-2, 6.31d-3, 1.27d-2, 1.0d-2, 1.0d-4, 3.63d-5, 2.22d-5]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp4_recommended, energy

; From GOESN-ENG-029, Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 4
   
   xenergy = [245., 255., 345., 355.]*1.0d0

   gf = [1.0d-4, 1.0d-2, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp4_composite, energy

; From GOESN-ENG-029, Table 3-29 and Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 9
   
   xenergy = [157.0, 200., 218., 241., 255., $
      274., 345., 355., 521.]*1.0d0

   gf = [1.22d-6, 2.37d-6, 5.39d-7, 6.60d-4, 1.0d-2, $
      7.47d-3, 1.0d-2, 1.0d-4, 1.95d-5]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
;   endif else if energy ge xenergy[2] and energy le xenergy[3] then begin
;      return, 1.0d-2
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp5_recommended, energy

; From GOESN-ENG-029, Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 4
   
   xenergy = [345., 355., 795., 805.]*1.0d0

   gf = [1.0d-4, 1.0d-2, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

function gf_mp5_composite, energy

; From GOESN-ENG-029, Table 3-29 and Table 3-36
; energy = proton energy in keV
; geometrical factor in cm^2 sr

   energy = double(energy)
   
   nenergy = 9
   
   xenergy = [200., 218., 241., 274., $
      345., 355., 521., 795., 805.]*1.0d0

   gf = [6.9d-7, 2.35d-7, 1.23d-6, 1.02d-6, $
      1.0d-4, 1.0d-2, 8.81d-3, 1.0d-2, 1.0d-4]
   
   if energy lt xenergy[0] or energy gt xenergy[nenergy-1] then begin
      return, 0.0
   endif else begin
      return, 10.0^interpol(alog10(gf), xenergy, energy)
   endelse


end

;------------------------------------------------------------------------------------------

pro magpd_gf, energy

   energy = findgen(401)*5.0d0 ; 0-2000 keV
   gf = dblarr(401, 5)
;   energy = findgen(1001)*2.0d0 ; 0-2000 keV
;   gf = dblarr(1001, 5)

;   BpStr = 'Measured'   
;   for i = 0, 400 do begin
;      gf[i,0] = gf_mp1_measured(energy[i])
;      gf[i,1] = gf_mp2_composite(energy[i])
;      gf[i,2] = gf_mp3_composite(energy[i])
;      gf[i,3] = gf_mp4_composite(energy[i])
;      gf[i,4] = gf_mp5_composite(energy[i])
;   endfor

   BpStr = 'Ideal'
   for i = 0, 400 do begin     
      gf[i,0] = gf_MP1_recommended(energy[i])
      gf[i,1] = gf_MP2_recommended(energy[i])
      gf[i,2] = gf_MP3_recommended(energy[i]) 
      gf[i,3] = gf_MP4_recommended(energy[i])
      gf[i,4] = gf_MP5_recommended(energy[i])
   endfor
        
   set_plot, 'win'

   colorcat12, color_scale
   !p.multi = [0,1,1]

   window, 10
   
   xrng = [60, 2000]
   yrng = [1e-7, 1e-1]
   
   plot, energy, gf[*,0], /nodata, /xlog, /ylog, xr = xrng, yr = yrng, xstyle = 1, ystyle = 1, $
      xtitle = 'Proton Energy (keV)', ytitle = 'Geometrical Factor (cm!U2!N sr)', $
      title = BpStr + ' MAGPD Geometrical Factors'

   for j = 0, 4 do oplot, energy, gf[*,j], color = color_scale[(2*j+1)], thick = 2

   Img = TVRD(0, TRUE = 1)
   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\bowtie\magpd_gf_loglog.png'
   WRITE_PNG, PlotnamePNG, Img
   
   window, 11
   
   xrng = [10, 2000]
   yrng = [0, 1.8e-2]
   
   plot, energy, gf[*,0], /nodata, /xlog, xr = xrng, yr = yrng, xstyle = 1, ystyle = 1, $
      xtitle = 'Proton Energy (keV)', ytitle = 'Geometrical Factor (cm!U2!N sr)', $
      title = BpStr + ' MAGPD Geometrical Factors'

   for j = 0, 4 do oplot, energy, gf[*,j], color = color_scale[(2*j+1)], thick = 2

   Img = TVRD(0, TRUE = 1)
   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\bowtie\magpd_gf_linlog_'+BpStr+'.png'
   WRITE_PNG, PlotnamePNG, Img
   
;   Img = TVRD(0, TRUE = 1)
;   PlotnamePNG = 'C:\Data\GOES-R\Phase 2\PAD inversion\weighting_function.png'
;   PlotnamePNG = 'C:\Data\GOES-R\PAD inversion\weighting_function.png'

;   WRITE_PNG, PlotnamePNG, Img

; stop
      
end

