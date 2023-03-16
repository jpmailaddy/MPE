;------------- functions specific to a relativistic Maxwellian ---------------
;-----------------------------------------------------------------------------

function relmaxwellp, q, energy

; Produces relativistic Maxwellian as function of q
; input energy in keV
   
   COMPILE_OPT IDL2
   
   q = double(q)
   
   ;rest mass
   mass_p_keV = 938272.013d0 
   
   flux = energy*(1.0d0 + energy/mass_p_keV/2.)*exp(q[0] + q[1]*energy)
   
   return, flux
   
end


;-----------------------------------------------------------------------------

function relmaxwellp_gradient, q, energy

; Produces first derivatives of relativistic Maxwellian with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   ;rest mass
   mass_p_keV = 938272.013d0 
   
   dfdq1 = energy*(1.0d0 + energy/mass_p_keV/2.)*exp(q[0] + q[1]*energy)
   dfdq2 = energy*dfdq1
   
   if n_elements(energy) eq 1 then begin
      return, [dfdq1, dfdq2]
   endif else begin   
      return, [transpose(dfdq1), transpose(dfdq2)]
   endelse
   
end

;-----------------------------------------------------------------------------

function relmaxwellp_curvature, q, energy

; Produces second derivatives of relativistic Maxwellian with respect to q
; input energy in keV
   
   COMPILE_OPT IDL2

   q = double(q)

   ;rest mass
   mass_p_keV = 938272.013d0 
   
   d2fdq12 = energy*(1.0d0 + energy/mass_p_keV/2.)*exp(q[0] + q[1]*energy)
   d2fdq1dq2 = energy*d2fdq12
   d2fdq22 = energy*d2fdq1dq2

   if n_elements(energy) eq 1 then begin
      return, [d2fdq12, d2fdq1dq2, d2fdq22]
   endif else begin   
      return, [transpose(d2fdq12), transpose(d2fdq1dq2), transpose(d2fdq22)]
   endelse
   
end

;----------------------------------------------------------------------------------

function q_first_guess_RMp, jflux, chan_energy

; first guess for a relativistic Maxwellian (RM)
; chan_energy is channel energy in keV  

   c_cps = 299792458.0d2 
   mass_p_keV = 938272.013d0 
   
   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP
   
   ERelSum = chan_energy*(chan_energy + 2.0d0*mass_p_keV) ; in kev^2

   ;PSD = c_cps*c_cps*jflux/ERelSum ; in (keV*s)^-3
   PSD = c_cps*c_cps*tflux/ERelSum ; in (keV*s)^-3
      
   FitPar = LINFIT( chan_energy, alog(PSD)) 

   q1 = alog(2.0*mass_p_keV*exp(FitPar[0])/c_cps/c_cps)

   Q0 = [q1, FitPar[1]]

   return, Q0
end

pro np_T_from_q, q, n, T

; Only possible for a Maxwellian.  This is for a non-relativistic Maxwellian
   COMPILE_OPT IDL2

   mass_p_keV = 938272.013d0 
   c_cps = 299792458.0d2
   
   modbess
   
   T = -1.0d0/q[1]
   
;   BetaRatio = -q[1]*mass_e_keV
;   
;   n = 2.0d0*!dpi*T^2*BetaRatio*exp(BetaRatio)*modbess_k2(BetaRatio)*exp(q[0])/c_cps
   
;   stop
   
end


;------------------ functions specific to a Double Relativistic Maxwellian ---------------------
;----------------------------------------------------------------------------------
function q_first_guess_DMp, jflux, chan_energy

; first guess for a relativistic Maxwellian (RM)
; chan_energy is channel energy in keV  
   c_cps = 299792458.0d2
   mass_p_keV = 938272.013d0

   tflux = jflux ;EDP, using tflux (temporary flux) so as not to change jflux. This function has issues with zero values.
   tflux[where(tflux eq 0.)] = 0.00390625 ;below instrument sensitivy with error (.0625^2.)  ;EDP

   ERelSum = chan_energy*(chan_energy + 2.0d0*mass_p_keV) ; in kev^2

   ;PSD = c_cps*c_cps*jflux/ERelSum ; in (keV*s)^-3
   PSD = c_cps*c_cps*tflux/ERelSum ; in (keV*s)^-3

   FitPar1 = LINFIT( chan_energy[0:2], alog(PSD[0:2]))
   q1 = alog(2.0*mass_p_keV*exp(FitPar1[0])/c_cps/c_cps)

   FitPar2 = LINFIT( chan_energy[2:4], alog(PSD[2:4]))
   q3 = alog(2.0*mass_p_keV*exp(FitPar2[0])/c_cps/c_cps)

   Q0 = [q1, FitPar1[1],q3,FitPar2[1]]
   return, Q0
end
;-----------------------------------------------------------------------------

pro drelmaxwellp, X, A, F, pder
;Produces double maxwellian as function of 4 free variables q.
;This procedure is for use with CURVEFIT.
;E = energy [keV]
;A = q = vector of free variables first guess.
;f = first guess of solution
;pder = vector of derivatives

   COMPILE_OPT IDL2

   E0 = 938272.013d0 

   E = X
   q = A

   base = E*(1.+(E/E0/2.))
   e1 = exp(q[0]+(q[1]*E))
   e2 = exp(q[2]+(q[3]*E))
   F = base*(e1+e2)

   ;If the procedure is called with four parameters, calculate the partial derivatives.

   if N_PARAMS() ge 4 then $
    pder = drelmaxwellp_gradiant(q,E)
    ;stop
end
;-----------------------------------------------------------------------------

function fdrelmaxwellp, q, energy 
;Produces double maxwellian as function of 4 free variables q.
;this is for use with call_function
;E = energy [keV]
;q = vector of free variables first guess.
;f = first guess of solution
;qder = vector of derivatives

   COMPILE_OPT IDL2

   E0 = 938272.013d0

   E = energy

   base = E*(1.+(E/E0/2.))
   e1 = exp(q[0]+(q[1]*E))
   e2 = exp(q[2]+(q[3]*E))
   F = base*(e1+e2)
   return, F
end
;-----------------------------------------------------------------------------

function drelmaxwellp_gradiant,q,energy
;for use with CURVEFIT to get pver
   q = double(q)
   E = double(energy)
   E0 = 938272.013d0 

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
     ;return, [transpose(dfdq1),transpose(dfdq2),transpose(dfdq3),transpose(dfdq4)]
     return, [[dfdq1],[dfdq2],[dfdq3],[dfdq4]]
   endelse
end
;-----------------------------------------------------------------------------

function fdrelmaxwellp_gradient,q,energy
;for use with call_function in calculating dl/dq
   q = double(q)
   E = double(energy)
   E0 = 938272.013d0

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

function fdrelmaxwellp_curvature, q, energy

; Produces second derivatives of relativistic Maxwellian with respect to q
; input energy in keV

   COMPILE_OPT IDL2

   q = double(q)
   E = double(energy)
   n = N_ELEMENTS(energy)
   E0 = 938272.013d0

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
     
   ;FitPar = LINFIT( alog(chan_energy), alog(jflux)) 
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
     
   ;Q0 = LINFIT( chan_energy, alog(jflux))
   Q0 = LINFIT( chan_energy, alog(tflux))
   
   return, Q0
end

;--------------------------- general functions -------------------------------
;-----------------------------------------------------------------------------

pro magpd_weighting_function, MEASURED = measured

; Calculates 5 row x Nenergy column weighting function matrix from MAGED geometrical factors
; Creates evenly-spaced, linear energy grid in keV.

; Called only once, at beginning of run-time.

; Inputs:
; dt - scalar in seconds

; Outputs:
; energy - Nenergy-element vector in keV
; Kwfp - (Nenergy x 5) weighting function matrix in cm^2 sr keV s
; Note: this weighting function gives counts (not count rate) when it is multiplied by differential flux

   COMPILE_OPT IDL2

   common weighting, dt, energy, Kwfe, Kwfp
   
   if ~keyword_set(measured) then measured = 0
      
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
   
; Define Kwfp
   Kwfp = fltarr(Nenergy, 5)
   
   if measured eq 0 then begin
   for i = 0, Nenergy-1 do begin
      Kwfp[i,0] = gf_mp1_recommended(energy[i])*dE[i]
      Kwfp[i,1] = gf_mp2_recommended(energy[i])*dE[i]
      Kwfp[i,2] = gf_mp3_recommended(energy[i])*dE[i]
      Kwfp[i,3] = gf_mp4_recommended(energy[i])*dE[i]         
      Kwfp[i,4] = gf_mp5_recommended(energy[i])*dE[i]
   endfor
   endif else begin
   for i = 0, Nenergy-1 do begin
      Kwfp[i,0] = gf_mp1_measured(energy[i])*dE[i]
      Kwfp[i,1] = gf_mp2_composite(energy[i])*dE[i]
      Kwfp[i,2] = gf_mp3_composite(energy[i])*dE[i]
      Kwfp[i,3] = gf_mp4_composite(energy[i])*dE[i]
      Kwfp[i,4] = gf_mp5_composite(energy[i])*dE[i]
   endfor
   endelse
      
   Kwfp *= dt
   
   print, 'MAGPD weighting function calculated....'
   
end
                           
;-----------------------------------------------------------------------------

function forward_model_magpd, q, jchan, spectral_function

; For a single channel, returns counts for an analytical function whose name
; is 'spectral_function' with parameter vector q

; Common inputs:
; energy - Nenergy-element vector in keV
; Kwfp - (Nenergy x 5) weighting function matrix in cm^2 sr keV s

   COMPILE_OPT IDL2

   common weighting, dt, energy, Kwfe, Kwfp
   
   Nenergy = n_elements(energy)
   ;stop 
   cts = 0.0d0
   for i = 0, Nenergy-1 do begin
      cts += call_function(spectral_function, q, energy[i])*Kwfp[i, jchan]
   endfor

   return, cts

end

;-----------------------------------------------------------------------------

function forward_model_magpd2, q, jchan, spectral_function

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
   cts = Kwfp##f

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
   
   ;if min(Varflux) lt 0 then stop
;   stop
   
   sd_logflux = sqrt(VarFlux); equation 42 in the O'Brien document.
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
   ;if min(Varflux) lt 0 then stop
   ;sometimes VarFlux is less than zero due to negative values from the covariance matrix. I am going to absolute value the varFlux. I think this is valid, but I should check with Juan.
   VarFlux = abs(VarFlux)
;   stop

   sd_logflux = sqrt(VarFlux)
   sd_flux = flux*sd_logflux ; convert std dev of ln(flux) to std dev of flux

end

;-----------------------------------------------------------------------------



function penalty_function_magpd_gradient, q, spectral_function, $
   spectral_function_1stderiv

; returns two-element vector containing dl/dq1 and dl/dq2
; 
; Kwfp (5r x 200c)
; energy (200 elements)
   common weighting, dt, energy, Kwfe, Kwfp

; ell (3r x 5c)   
   ell = penalty_function_magpd(q, spectral_function)

; dfdq (200r x 2c)   
   dfdq = call_function(spectral_function_1stderiv, q, energy) 
   
; FwdDeriv: forward model derivatives for the 5 channels: 5r x 2c
   FwdDeriv = Kwfp##dfdq

; dldq: 5r x 2c  
   dldq_1 = ell[0:4,1]*FwdDeriv[0,*]
   dldq_2 = ell[0:4,1]*FwdDeriv[1,*]
   if (size(FwdDeriv))[1] eq 4 then begin
     dldq_3 = ell[0:4,1]*FwdDeriv[2,*]
     dldq_4 = ell[0:4,1]*FwdDeriv[3,*]
     ;stop
     return, [total(dldq_1), total(dldq_2), total(dldq_3), total(dldq_4)]
   endif else return, [total(dldq_1), total(dldq_2)]
   
end

;-----------------------------------------------------------------------------

function penalty_function_magpd_hessian, q, spectral_function, $
   spectral_function_1stderiv, spectral_function_2ndderiv

; returns four-element array containing d2l/dq2
; Kwfp (5r x 200c)
; energy (200 elements)
   common weighting, dt, energy, Kwfe, Kwfp

; ell (3r x 5c)    
   ell = penalty_function_magpd(q, spectral_function)

; dfdq (200r x 2c)   
   dfdq = call_function(spectral_function_1stderiv, q, energy)

; d2fdq2 (200r x 3c)     
   d2fdq2 = call_function(spectral_function_2ndderiv, q, energy)

; FwdDeriv: forward model derivatives for the 5 channels: 5r x 2c
   FwdDeriv = Kwfp##dfdq
   d1size = size(FwdDeriv)
; Fwd2ndDeriv: forward model 2nd derivatives for the 5 channels: 5r x 3c
   Fwd2ndDeriv = Kwfp##d2fdq2

; d2ldq2_summand1: 5r x 3c  
   if d1size[1] eq 4 then d2ldq2_summand1 = dblarr(10, 5) else d2ldq2_summand1 = dblarr(3, 5)
   d2ldq2_summand1[0, *] = ell[0:4,2]*FwdDeriv[0,*]^2   ; q1*q1
   d2ldq2_summand1[1, *] = ell[0:4,2]*FwdDeriv[0,*]*FwdDeriv[1,*] ;q1*q2
   d2ldq2_summand1[2, *] = ell[0:4,2]*FwdDeriv[1,*]^2  ;q2*q2
   ;need to add in more for q3 and q4 in evend of double maxwellian spectra
   if d1size[1] eq 4 then begin
     d2ldq2_summand1[3, *] = ell[0:4,2]*FwdDeriv[0,*]*FwdDeriv[2,*] ;q1*q3
     d2ldq2_summand1[4, *] = ell[0:4,2]*FwdDeriv[0,*]*FwdDeriv[3,*] ;q1*q4
     d2ldq2_summand1[5, *] = ell[0:4,2]*FwdDeriv[1,*]*FwdDeriv[2,*] ;q2*q3
     d2ldq2_summand1[6, *] = ell[0:4,2]*FwdDeriv[1,*]*FwdDeriv[3,*] ;q2*q4
     d2ldq2_summand1[7, *] = ell[0:4,2]*FwdDeriv[2,*]^2   ; q3*q3
     d2ldq2_summand1[8, *] = ell[0:4,2]*FwdDeriv[2,*]*FwdDeriv[3,*] ;q3*q4
     d2ldq2_summand1[9, *] = ell[0:4,2]*FwdDeriv[3,*]^2   ; q4*q4
   endif 


; d2ldq2_summand2: 5r x 3c  
   if d1size[1] eq 4 then d2ldq2_summand2 = dblarr(10,5) else d2ldq2_summand2 = dblarr(3, 5)
   d2ldq2_summand2[0, *] = ell[0:4,1]*Fwd2ndDeriv[0,*] ;dq1 dq1      
   d2ldq2_summand2[1, *] = ell[0:4,1]*Fwd2ndDeriv[1,*] ;dq1 dq2
   d2ldq2_summand2[2, *] = ell[0:4,1]*Fwd2ndDeriv[2,*] ;dq2 dq2
   if d1size[1] eq 4 then begin
     d2ldq2_summand2[3, *] = ell[0:4,1]*Fwd2ndDeriv[3,*] ;dq1 dq3
     d2ldq2_summand2[4, *] = ell[0:4,1]*Fwd2ndDeriv[4,*] ;dq1 dq4
     d2ldq2_summand2[5, *] = ell[0:4,1]*Fwd2ndDeriv[5,*] ;dq2 dq3
     d2ldq2_summand2[6, *] = ell[0:4,1]*Fwd2ndDeriv[6,*] ;dq2 dq4
     d2ldq2_summand2[7, *] = ell[0:4,1]*Fwd2ndDeriv[7,*] ;dq3 dq3
     d2ldq2_summand2[8, *] = ell[0:4,1]*Fwd2ndDeriv[8,*] ;dq3 dq4
     d2ldq2_summand2[9, *] = ell[0:4,1]*Fwd2ndDeriv[9,*] ;dq4 dq4
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

function penalty_function_magpd, q, spectral_function
;this needs to change to be more dynamic depending on number of terms.
; Output: ell, 6 x 3 element array
; first row: individual ell's at each y, and sum (6 terms)
; second row: individual first derivatives (dell/dlambda), last element is not defined
; third row: individual second derivatives (d2ell/dlambda2), last element is not defined
; (There are no cross terms in the second derivatives)

   common counts, y, dy
  
   Nchannels = 5 ;EDP
 
   ell = fltarr(Nchannels+1,3) ;adapted for number of channels. for protons P1-P5 and sum will equal 6. first dimension is number of channels. second dimension is l, dl/dlambda, and d2l,dlambda2.
   
   lambda_vec = forward_model_magpd2(q, i, spectral_function)
   for i = 0, 4 do begin      
       ;lambda1 = forward_model_magpd(q, i, spectral_function) ;test ;EDP
       lambda = lambda_vec[i]
       ;stop
      if y[i] lt dy^(-2) then begin ; Poisson (counting statistics)
         
         ell[i,0] = lambda - y[i]*alog(lambda)
         ell[i,1] = 1.0 - y[i]/lambda
         ell[i,2] = y[i]/lambda^2
      
      endif else begin ; Gaussian (calibration)
      
         ell[i,0] = 0.5d0*((alog(y[i])-alog(lambda))/dy)^2
         ell[i,1] = (alog(lambda)-alog(y[i]))/(dy^2)/lambda
         ell[i,2] = (1.0 + alog(y[i])-alog(lambda))/(lambda*dy)^2
      
      endelse
      
;      print, y[i], lambda, ell[i,0]
   
   endfor
   ell[5,0] = total(ell[0:4,0])

   return, ell
   
end




