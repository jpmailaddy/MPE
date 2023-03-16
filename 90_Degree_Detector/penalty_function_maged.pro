;------------- functions specific to a relativistic Maxwellian ---------------
;-----------------------------------------------------------------------------

function relmaxwell, q, energy

; Produces relativistic Maxwellian as function of q
; input energy in keV
   
   COMPILE_OPT IDL2
   
   q = double(q)
   
   ;rest mass
   mass_e_keV = 510.998910d0
   
   flux = energy*(1.0d0 + energy/mass_e_keV/2.)*exp(q[0] + q[1]*energy)
   
   return, flux
   
end


;-----------------------------------------------------------------------------

function relmaxwell_gradient, q, energy

; Produces first derivatives of relativistic Maxwellian with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   ;rest mass
   mass_e_keV = 510.998910d0
   
   dfdq1 = energy*(1.0d0 + energy/mass_e_keV/2.)*exp(q[0] + q[1]*energy)
   dfdq2 = energy*dfdq1
   
   if n_elements(energy) eq 1 then begin
      return, [dfdq1, dfdq2]
   endif else begin   
      return, [transpose(dfdq1), transpose(dfdq2)]
   endelse
   
end

;-----------------------------------------------------------------------------

function relmaxwell_curvature, q, energy

; Produces second derivatives of relativistic Maxwellian with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   ;rest mass
   mass_e_keV = 510.998910d0
   
   d2fdq12 = energy*(1.0d0 + energy/mass_e_keV/2.)*exp(q[0] + q[1]*energy)
   d2fdq1dq2 = energy*d2fdq12
   d2fdq22 = energy*d2fdq1dq2

   if n_elements(energy) eq 1 then begin
      return, [d2fdq12, d2fdq1dq2, d2fdq22]
   endif else begin   
      return, [transpose(d2fdq12), transpose(d2fdq1dq2), transpose(d2fdq22)]
   endelse
   
end

;----------------------------------------------------------------------------------

function q_first_guess_RM, jflux, chan_energy

; first guess for a relativistic Maxwellian (RM)
; chan_energy is channel energy in keV  

   c_cps = 299792458.0d2
   keV_to_erg = 1.602176487d-9   
   mass_e_keV = 510.998910d0
     
   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP
 
   ERelSum = chan_energy*(chan_energy + 2.0d0*mass_e_keV) ; in kev^2

   PSD = c_cps*c_cps*tflux/ERelSum ; in (keV*s)^-3
      
   FitPar = LINFIT( chan_energy, alog(PSD)) 

   q1 = alog(2.0*mass_e_keV*exp(FitPar[0])/c_cps/c_cps)

   Q0 = [q1, FitPar[1]]
   return, Q0
end

;----------------------------------------------------------------------------------

pro n_T_from_q, q, n, T

; Only possible for a Maxwellian.  This is for a relativistic Maxwellian
   COMPILE_OPT IDL2

   mass_e_keV = 510.998910d0 
   c_cps = 299792458.0d2
   
   modbess
   
   T = -1.0d0/q[1]
   
   BetaRatio = -q[1]*mass_e_keV
   
   n = 2.0d0*!dpi*T^2*BetaRatio*exp(BetaRatio)*modbess_k2(BetaRatio)*exp(q[0])/c_cps
   
;   stop
   
end

;------------------ functions specific to a double maxwellian -------------
;----------------------------------------------------------------------------------

function q_first_guess_DM, jflux, chan_energy

; first guess for a relativistic Double Maxwellian (DM)
; chan_energy is channel energy in keV  

   c_cps = 299792458.0d2
   mass_e_keV = 510.998910d0

   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP

   ERelSum = chan_energy*(chan_energy + 2.0d0*mass_e_keV) ; in kev^2

   PSD = c_cps*c_cps*tflux/ERelSum ; in (keV*s)^-3

   FitPar1 = LINFIT( chan_energy[0:2], alog(PSD[0:2]))
   q1 = alog(2.0*mass_e_keV*exp(FitPar1[0])/c_cps/c_cps)

   FitPar2 = LINFIT( chan_energy[1:3], alog(PSD[1:3]))
   q3 = alog(2.0*mass_e_keV*exp(FitPar2[0])/c_cps/c_cps)

   Q0 = [q1, FitPar1[1], q3, FitPar2[1]]
   return, Q0
end
;-----------------------------------------------------------------------------
function fdrelmaxwell, q, energy

   COMPILE_OPT IDL2

   q = double(q)
   E = double(energy)
   E0 = 510.998910d0

   base = E*(1.+(E/E0/2.))
   e1 = exp(q[0]+(q[1]*E))
   e2 = exp(q[2]+(q[3]*E))
   F = base*(e1+e2)
   return, F
end
;-----------------------------------------------------------------------------

function fdrelmaxwell_gradient,q,energy
   q = double(q)
   E = double(energy)
   E0 = 510.998910d0
  
   base = E*(1.+(E/E0/2.))
   e1 = exp(q[0]+(q[1]*E))
   e2 = exp(q[2]+(q[3]*E))

   dfdq1 = base*e1
   dfdq2 = E*dfdq1
   dfdq3 = base*e2
   dfdq4 = E*dfdq3

   if N_ELEMENTS(E) eq 1 then begin
     return, [dfdq1, dfdq2, dfdq3, dfdq4]
   endif else begin
     return, [transpose(dfdq1),transpose(dfdq2),transpose(dfdq3),transpose(dfdq4)]
   endelse
end
;-----------------------------------------------------------------------------

function fdrelmaxwell_curvature,q,energy
   COMPILE_OPT IDL2

   q = double(q)
   E = double(energy)
   n = N_ELEMENTS(energy)
   E0 = 510.998910d0
   
   base = E*(1.+(E/E0/2.))
   e1 = exp(q[0]+(q[1]*E))
   e2 = exp(q[2]+(q[3]*E))
   
   d2fdq12 = base*e1
   d2fdq1dq2 = E*d2fdq12
   d2fdq22 = E*d2fdq1dq2
   d2fdq1dq3 = dblarr(n)
   d2fdq1dq4 = dblarr(n)
   d2fdq2dq3 = dblarr(n)
   d2fdq2dq4 = dblarr(n) 
   d2fdq32 = base*e2 
   d2fdq3dq4 = E*d2fdq32
   d2fdq42 = E*d2fdq3dq4

   if n_elements(energy) eq 1 then begin
      return, [d2fdq12, d2fdq1dq2, d2fdq22, d2fdq1dq3, d2fdq1dq4, d2fdq2dq3, d2fdq2dq4, d2fdq32, d2fdq3dq4, d2fdq42]
   endif else begin   
      return, [transpose(d2fdq12), transpose(d2fdq1dq2), transpose(d2fdq22),transpose(d2fdq1dq3), transpose(d2fdq1dq4), transpose(d2fdq2dq3), transpose(d2fdq2dq4), transpose(d2fdq32), transpose(d2fdq3dq4), transpose(d2fdq42)]
   endelse
end

;------------------ functions specific to a power law ---------------------
;-----------------------------------------------------------------------------

function powerlaw, q, energy

; Produces power law as function of q
; input energy in keV
   
   COMPILE_OPT IDL2
   
   q = double(q)
   
   flux = exp(q[0] - q[1]*alog(energy))
   
   return, flux
   
end

;-----------------------------------------------------------------------------

function powerlaw_gradient, q, energy

; Produces first derivatives of power law with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)
   
   dfdq1 = exp(q[0] - q[1]*alog(energy))
   dfdq2 = -alog(energy)*dfdq1
   
   if n_elements(energy) eq 1 then begin
      return, [dfdq1, dfdq2]
   endif else begin   
      return, [transpose(dfdq1), transpose(dfdq2)]
   endelse
   
end

;-----------------------------------------------------------------------------

function powerlaw_curvature, q, energy

; Produces second derivatives of power law with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   d2fdq12   = exp(q[0] - q[1]*alog(energy))
   d2fdq1dq2 = -alog(energy)*d2fdq12
   d2fdq22   = -alog(energy)*d2fdq1dq2

   if n_elements(energy) eq 1 then begin
      return, [d2fdq12, d2fdq1dq2, d2fdq22]
   endif else begin   
      return, [transpose(d2fdq12), transpose(d2fdq1dq2), transpose(d2fdq22)]
   endelse
   
end

;----------------------------------------------------------------------------------

function q_first_guess_PL, jflux, chan_energy

   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP


; chan_energy is channel energy in keV  
     
   FitPar = LINFIT( alog(chan_energy), alog(tflux)) 

   Q0 = [FitPar[0], -FitPar[1]]
;   print, q0
   
;   stop
   
   return, Q0
end

;----------- functions specific to an energy exponential  --------------------
;-----------------------------------------------------------------------------

function energyexponential, q, energy

; Produces energy exponential as function of q
; input energy in keV
   
   COMPILE_OPT IDL2
   
   q = double(q)
   
   flux = exp(q[0] + q[1]*energy)
   
   return, flux
   
end

;-----------------------------------------------------------------------------

function energyexponential_gradient, q, energy

; Produces first derivatives of energy exponential with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)
   
   dfdq1 = exp(q[0] + q[1]*energy)
   dfdq2 = energy*dfdq1
   
   if n_elements(energy) eq 1 then begin
      return, [dfdq1, dfdq2]
   endif else begin   
      return, [transpose(dfdq1), transpose(dfdq2)]
   endelse
   
end

;-----------------------------------------------------------------------------

function energyexponential_curvature, q, energy

; Produces second derivatives of energy exponential with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   d2fdq12   = exp(q[0] + q[1]*energy)
   d2fdq1dq2 = energy*d2fdq12
   d2fdq22   = energy*d2fdq1dq2

   if n_elements(energy) eq 1 then begin
      return, [d2fdq12, d2fdq1dq2, d2fdq22]
   endif else begin   
      return, [transpose(d2fdq12), transpose(d2fdq1dq2), transpose(d2fdq22)]
   endelse
   
end

;----------------------------------------------------------------------------------

function q_first_guess_EE, jflux, chan_energy

   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP

; chan_energy is channel energy in keV  
     
   Q0 = LINFIT( chan_energy, alog(tflux))
   
   return, Q0
end

;--------------------------- general functions -------------------------------
;-----------------------------------------------------------------------------

pro maged_weighting_function

; Calculates 5 row x Nenergy column weighting function matrix from MAGED geometrical factors
; Creates evenly-spaced, linear energy grid in keV.

; Called only once, at beginning of run-time.

; Inputs:
; dt - scalar in seconds

; Outputs:
; energy - Nenergy-element vector in keV
; Kwfe - (Nenergy x 5) weighting function matrix in cm^2 sr keV s
; Note: this weighting function gives counts (not count rate) when it is multiplied by differential flux

   COMPILE_OPT IDL2

   common weighting, dt, energy, Kwfe, Kwfp
      
; Define the energy scale
   Nenergy = 200
   DeltaE = float(2000/Nenergy)
   energy = findgen(Nenergy)*DeltaE + DeltaE ; 0-2000 keV
   
; Define dE for trapezoidal integration following O'Brien (2010) consistent with Twomey (1977)
   dE = fltarr(Nenergy)
   dE[0] = (energy[1]-energy[0])
   dE[Nenergy-1] = (energy[Nenergy-1]-energy[Nenergy-2])
   for i = 1, Nenergy-2 do dE[i] = (energy[i+1]-energy[i-1])
   dE *= 0.5d0
   
; Define Kwfe
   Kwfe = fltarr(Nenergy, 5)
   for i = 0, Nenergy-1 do begin
      Kwfe[i,0] = gf_me1_recommended(energy[i])*dE[i]
      Kwfe[i,1] = gf_me2_recommended(energy[i])*dE[i]
      Kwfe[i,2] = gf_me3_recommended(energy[i])*dE[i]
      Kwfe[i,3] = gf_me4_recommended(energy[i])*dE[i]         
      Kwfe[i,4] = gf_me5_recommended(energy[i])*dE[i]
   endfor
   
   Kwfe *= dt
   
   print, 'MAGED weighting function calculated....'
   
end
                           
;-----------------------------------------------------------------------------

function forward_model_maged, q, jchan, spectral_function, test

; For a single channel, returns counts for an analytical function whose name
; is 'spectral_function' with parameter vector q

; Common inputs:
; energy - Nenergy-element vector in keV
; Kwfe - (Nenergy x 5) weighting function matrix in cm^2 sr keV s

   COMPILE_OPT IDL2

   common weighting, dt, energy, Kwfe, Kwfp
   
   Nenergy = n_elements(energy)
   ;if test eq 1 then stop
   cts = 0.0d0
   for i = 0, Nenergy-1 do begin
      cts += call_function(spectral_function, q, energy[i])*Kwfe[i, jchan]
   endfor
   ;if test eq 1 then stop
   return, cts

end

;-----------------------------------------------------------------------------

function forward_model_maged2, q, spectral_function;, test

; For a single channel, returns counts for an analytical function whose name
; is 'spectral_function' with parameter vector q

; Common inputs:
; energy - Nenergy-element vector in keV
; Kwfp - (Nenergy x 5) weighting function matrix in cm^2 sr keV s

   COMPILE_OPT IDL2

   common weighting, dt, energy, Kwfe, Kwfp

   Nenergy = n_elements(energy)

   ;cts = 0.0d0
   f = call_function(spectral_function,q,energy)
   cts = Kwfe##f
   ;if test eq 1 then stop
   return, cts
   ;for i = 0, Nenergy-1 do begin
   ;   cts += call_function(spectral_function, q, energy[i])*Kwfp[i, jchan]
   ;endfor

   ;return, cts

end


;-----------------------------------------------------------------------------
pro logflux_standard_deviation_2, q, covqq, energy, spectral_function, $
   spectral_function_1stderiv, sd_logflux, sd_flux

; For two-parameter fit.
; Q contains the two retrieved fit parameters q1 and q2 for the relativistic maxwellian
; d2ldq2 is the 2 x 2 covariance matrix of q1 and q2
; input energy in keV

   flux = call_function(spectral_function, q, energy)
; dfdq (200r x 2c)
   dfdq = call_function(spectral_function_1stderiv, q, energy)
   
   dfdq1_over_flux = dfdq[0,*]/flux
   dfdq2_over_flux = dfdq[1,*]/flux
   
   Var11 = dfdq1_over_flux*covqq[0,0]*dfdq1_over_flux
   Var12 = dfdq1_over_flux*covqq[0,1]*dfdq2_over_flux
   Var21 = dfdq2_over_flux*covqq[1,0]*dfdq1_over_flux
   Var22 = dfdq2_over_flux*covqq[1,1]*dfdq2_over_flux   
   
   VarFlux = Var11 + Var12 + Var21 + Var22 ; really, variance of ln of flux
   
;   stop
   
   sd_logflux = sqrt(VarFlux)
   sd_flux = flux*sd_logflux ; convert std dev of ln(flux) to std dev of flux
   
end

;-----------------------------------------------------------------------------

pro logflux_standard_deviation_4, q, covqq, energy, spectral_function, $
   spectral_function_1stderiv, sd_logflux, sd_flux

; For four-parameter fit.
; Q contains the four retrieved fit parameters q1,q2, q3, and q4 for the double relativistic maxwellian
; d2ldq2 is the 4 x 4 covariance matrix of q1 and q2
; input energy in keV

   flux = call_function(spectral_function, q, energy)
; dfdq (200r x 2c)
   dfdq = call_function(spectral_function_1stderiv, q, energy)

   dfdq1_over_flux = dfdq[0,*]/flux
   dfdq2_over_flux = dfdq[1,*]/flux
   dfdq3_over_flux = dfdq[2,*]/flux
   dfdq4_over_flux = dfdq[3,*]/flux

   Var11 = dfdq1_over_flux*covqq[0,0]*dfdq1_over_flux
   Var12 = dfdq1_over_flux*covqq[0,1]*dfdq2_over_flux
   Var13 = dfdq1_over_flux*covqq[0,2]*dfdq3_over_flux
   Var14 = dfdq1_over_flux*covqq[0,3]*dfdq4_over_flux
   Var21 = dfdq2_over_flux*covqq[1,0]*dfdq1_over_flux
   Var22 = dfdq2_over_flux*covqq[1,1]*dfdq2_over_flux
   Var23 = dfdq2_over_flux*covqq[1,2]*dfdq3_over_flux
   Var24 = dfdq2_over_flux*covqq[1,3]*dfdq4_over_flux
   Var31 = dfdq1_over_flux*covqq[2,0]*dfdq1_over_flux
   Var32 = dfdq1_over_flux*covqq[2,1]*dfdq2_over_flux
   Var33 = dfdq1_over_flux*covqq[2,2]*dfdq3_over_flux
   Var34 = dfdq1_over_flux*covqq[2,3]*dfdq4_over_flux
   Var41 = dfdq1_over_flux*covqq[3,0]*dfdq1_over_flux
   Var42 = dfdq1_over_flux*covqq[3,1]*dfdq2_over_flux
   Var43 = dfdq1_over_flux*covqq[3,2]*dfdq3_over_flux
   Var44 = dfdq1_over_flux*covqq[3,3]*dfdq4_over_flux

   VarFlux = Var11 + Var12 + Var13 + Var14 + Var21 + Var22 + Var23 + Var24 + Var31 + Var32 + Var33 + Var34 + Var41 + Var42 + Var43 + Var44 ; really, variance of ln of flux

;   stop
   VarFlux = abs(VarFlux) ;check with Juan to make sure this is Kosher.
   sd_logflux = sqrt(VarFlux)
   sd_flux = flux*sd_logflux ; convert std dev of ln(flux) to std dev of flux

end

;-----------------------------------------------------------------------------

function penalty_function_maged_gradient, q, spectral_function, $
   spectral_function_1stderiv

; returns two-element vector containing dl/dq1 and dl/dq2
; 
; Kwfe (5r x 200c)
; energy (200 elements)
   common weighting, dt, energy, Kwfe

; ell (3r x 5c)   
   ell = penalty_function_maged(q, spectral_function)

; dfdq (200r x 2c)   
   dfdq = call_function(spectral_function_1stderiv, q, energy) 
   
; FwdDeriv: forward model derivatives for the 5 channels: 5r x 2c
   FwdDeriv = Kwfe##dfdq

; dldq: 5r x 2c  
   dldq_1 = ell[0:3,1]*FwdDeriv[0,*]
   dldq_2 = ell[0:3,1]*FwdDeriv[1,*]  
   if (size(FwdDeriv))[1] eq 4 then begin
     dldq_3 = ell[0:3,1]*FwdDeriv[2,*]
     dldq_4 = ell[0:3,1]*FwdDeriv[3,*]
     return, [total(dldq_1), total(dldq_2), total(dldq_3), total(dldq_4)]
   endif else return, [total(dldq_1), total(dldq_2)]
   
end

;-----------------------------------------------------------------------------

function penalty_function_maged_hessian, q, spectral_function, $
   spectral_function_1stderiv, spectral_function_2ndderiv

; returns four-element array containing d2l/dq2
; Kwfe (5r x 200c)
; energy (200 elements)
   common weighting, dt, energy, Kwfe, Kwfp

; ell (3r x 5c)    
   ell = penalty_function_maged(q, spectral_function)

; dfdq (200r x 2c)   
   dfdq = call_function(spectral_function_1stderiv, q, energy)

; d2fdq2 (200r x 3c)     
   d2fdq2 = call_function(spectral_function_2ndderiv, q, energy)

; FwdDeriv: forward model derivatives for the 5 channels: 5r x 2c
   FwdDeriv = Kwfe##dfdq
   d1size = size(FwdDeriv)
   
; Fwd2ndDeriv: forward model 2nd derivatives for the 5 channels: 5r x 3c
   Fwd2ndDeriv = Kwfe##d2fdq2

; d2ldq2_summand1: 4r x 3c or 10c 
   if d1size[1] eq 4 then d2ldq2_summand1 = dblarr(10, 4) else d2ldq2_summand1 = dblarr(3, 4)
   d2ldq2_summand1[0, *] = ell[0:3,2]*FwdDeriv[0,*]^2   
   d2ldq2_summand1[1, *] = ell[0:3,2]*FwdDeriv[0,*]*FwdDeriv[1,*]
   d2ldq2_summand1[2, *] = ell[0:3,2]*FwdDeriv[1,*]^2 
   if d1size[1] eq 4 then begin
     d2ldq2_summand1[3, *] = ell[0:3,2]*FwdDeriv[0,*]*FwdDeriv[2,*] ;q1*q3
     d2ldq2_summand1[4, *] = ell[0:3,2]*FwdDeriv[0,*]*FwdDeriv[3,*] ;q1*q4
     d2ldq2_summand1[5, *] = ell[0:3,2]*FwdDeriv[1,*]*FwdDeriv[2,*] ;q2*q3
     d2ldq2_summand1[6, *] = ell[0:3,2]*FwdDeriv[1,*]*FwdDeriv[3,*] ;q2*q4
     d2ldq2_summand1[7, *] = ell[0:3,2]*FwdDeriv[2,*]^2   ; q3*q3
     d2ldq2_summand1[8, *] = ell[0:3,2]*FwdDeriv[2,*]*FwdDeriv[3,*] ;q3*q4
     d2ldq2_summand1[9, *] = ell[0:3,2]*FwdDeriv[3,*]^2   ; q4*q4
   endif

; d2ldq2_summand2: 4r x 3c or 10cforwar
   if d1size[1] eq 4 then d2ldq2_summand2 = dblarr(10,4) else d2ldq2_summand2 = dblarr(3, 4)
   d2ldq2_summand2[0, *] = ell[0:3,1]*Fwd2ndDeriv[0,*]      
   d2ldq2_summand2[1, *] = ell[0:3,1]*Fwd2ndDeriv[1,*]
   d2ldq2_summand2[2, *] = ell[0:3,1]*Fwd2ndDeriv[2,*]
   if d1size[1] eq 4 then begin
     d2ldq2_summand2[3, *] = ell[0:3,1]*Fwd2ndDeriv[3,*] ;dq1 dq3
     d2ldq2_summand2[4, *] = ell[0:3,1]*Fwd2ndDeriv[4,*] ;dq1 dq4
     d2ldq2_summand2[5, *] = ell[0:3,1]*Fwd2ndDeriv[5,*] ;dq2 dq3
     d2ldq2_summand2[6, *] = ell[0:3,1]*Fwd2ndDeriv[6,*] ;dq2 dq4
     d2ldq2_summand2[7, *] = ell[0:3,1]*Fwd2ndDeriv[7,*] ;dq3 dq3
     d2ldq2_summand2[8, *] = ell[0:3,1]*Fwd2ndDeriv[8,*] ;dq3 dq4
     d2ldq2_summand2[9, *] = ell[0:3,1]*Fwd2ndDeriv[9,*] ;dq4 dq4
   endif
 
   d2ldq2 = d2ldq2_summand1 + d2ldq2_summand2   
   
;   stop
  if d1size[1] eq 4 then begin
    return, [[total(d2ldq2[0,*]), total(d2ldq2[1,*]), total(d2ldq2[3,*]), total(d2ldq2[4,*])],$ ;q1
             [total(d2ldq2[1,*]), total(d2ldq2[2,*]), total(d2ldq2[5,*]), total(d2ldq2[6,*])],$ ;q2
             [total(d2ldq2[3,*]), total(d2ldq2[5,*]), total(d2ldq2[7,*]), total(d2ldq2[8,*])],$ ;q3
             [total(d2ldq2[4,*]), total(d2ldq2[6,*]), total(d2ldq2[8,*]), total(d2ldq2[9,*])]]  ;q4
  endif else begin
    return, [[total(d2ldq2[0,*]), total(d2ldq2[1,*])],[total(d2ldq2[1,*]), total(d2ldq2[2,*])]]
  endelse
end

;-----------------------------------------------------------------------------

function covariance_qq, d2ldq2

   cov = LA_INVERT( d2ldq2, /DOUBLE, status = InversionFlag)
   
   if InversionFlag gt 0 then begin
	print, 'd2l/dq2 is singular and the inverse could not be computed. '
	;stop
   endif
   
   return, cov
end

;-----------------------------------------------------------------------------

function penalty_function_maged, q, spectral_function;, test

; Output: ell, 5 x 3 element array
; first row: individual ell's at each y, and sum (5 terms)
; second row: individual first derivatives (dell/dlambda), last element is not defined
; third row: individual second derivatives (d2ell/dlambda2), last element is not defined
; (There are no cross terms in the second derivatives)

   common counts, y, dy
  
   Nchannels = 4 ;EDP ;E1,E2,E3,E4
   ;Nchannels = 3 ;EDP ;E1,E2,E3 or E4
 
   ell = fltarr(Nchannels+1,3)
   lambda_vec = forward_model_maged2(q, spectral_function);, test)
   ;if test eq 1 then stop
   for i = 0, Nchannels-1 do begin      
      lambda = forward_model_maged(q, i, spectral_function, test)
      ;lambda = lambda_vec[i]
      ;if test eq 1 then stop
      if y[i] lt dy^(-2) then begin ; Poisson (counting statistics)
         
         ell[i,0] = lambda - y[i]*alog(lambda)
         ell[i,1] = 1.0 - y[i]/lambda
         ell[i,2] = y[i]/lambda^2
      
      endif else begin ; Gaussian (calibration)
      
         ell[i,0] = 0.5d0*((alog(y[i])-alog(lambda))/dy)^2
         ell[i,1] = (alog(lambda)-alog(y[i]))/(dy^2)/lambda
         ell[i,2] = (1.0 + alog(y[i])-alog(lambda))/(lambda*dy)^2
      
      endelse
      
   endfor
   ell[Nchannels,0] = total(ell[0:Nchannels-1,0])
   ;if test eq 1 then stop
   return, ell
   
end




