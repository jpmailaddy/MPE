
; --------------------------------------------------------------------------------------------------------
; The following code runs the scripts, functions and procedures needed to both remove proton 
; contamination from the MEPED instrument electrons channels and calculate the differential flux spectrum
; from the POES satellites. 
; 
;                      Created by: Ethan D. Peck
;                      Modified by: Josh Pettit 
;                      Date: 7/19/2018
;                      
;                      
;                      Inputs: satellite and date
; 
; --------------------------------------------------------------------------------------------------------

@D:\POES_Data\MPE_Software\90_Degree_Detector\spectral_fits_90.pro

;stop

spectral_fits_jmp_v2_90

;ammendNCfiles_BLC_90
;ammendNCfiles_magField_90

end
