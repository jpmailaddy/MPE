;-----------------------------------------------------------------------------

function penalty_maged_exponential, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_maged(q, 'energyexponential')
   
   return, ell_vector[4]

end

;-----------------------------------------------------------------------------

function penalty_gradient_maged_exponential, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_maged_gradient(q, 'energyexponential', 'energyexponential_gradient') 
   
   return, dldq

end

;-----------------------------------------------------------------------------

function penalty_maged_powerlaw, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_maged(q, 'powerlaw')
   
   return, ell_vector[4]

end

;-----------------------------------------------------------------------------

function penalty_gradient_maged_powerlaw, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_maged_gradient(q, 'powerlaw', 'powerlaw_gradient') 
   
   return, dldq

end


;-----------------------------------------------------------------------------

function penalty_maged_relmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_maged(q, 'relmaxwell');,0)
   
   return, ell_vector[4]

end

;-----------------------------------------------------------------------------

function penalty_gradient_maged_relmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single MAGED telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_maged_gradient(q, 'relmaxwell', 'relmaxwell_gradient') 
   
   return, dldq

end

function penalty_maged_drelmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2

   common counts, y, dy

   y = double(y)
   dy = double(dy)
;       stop
   ell_vector = penalty_function_maged(q, 'fdrelmaxwell')

   return, ell_vector[4]

end

;-----------------------------------------------------------------------------

function penalty_gradient_maged_drelmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2

   common counts, y, dy

   y = double(y)
   dy = double(dy)

   dldq = penalty_function_maged_gradient(q, 'fdrelmaxwell', 'fdrelmaxwell_gradient')

   return, dldq

end



;-----------------------------------------------------------------------------

function penalty_magpd_exponential, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_magpd(q, 'energyexponential')
   
   return, ell_vector[5]

end

;-----------------------------------------------------------------------------

function penalty_gradient_magpd_exponential, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_magpd_gradient(q, 'energyexponential', 'energyexponential_gradient') 
   
   return, dldq

end

;-----------------------------------------------------------------------------

function penalty_magpd_powerlaw, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_magpd(q, 'powerlaw')
   
   return, ell_vector[5]

end

;-----------------------------------------------------------------------------

function penalty_gradient_magpd_powerlaw, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_magpd_gradient(q, 'powerlaw', 'powerlaw_gradient') 
   
   return, dldq

end


;-----------------------------------------------------------------------------

function penalty_magpd_relmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   ell_vector = penalty_function_magpd(q, 'relmaxwellp')
   
   return, ell_vector[5]

end

;-----------------------------------------------------------------------------

function penalty_gradient_magpd_relmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy
   
   y = double(y)
   dy = double(dy)

   dldq = penalty_function_magpd_gradient(q, 'relmaxwellp', 'relmaxwellp_gradient') 
   
   return, dldq

end

;-----------------------------------------------------------------------------

function penalty_magpd_drelmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by AMOEBA and DFPMIN.
; Derivatives (for DFPMIN) are calculated by a separate function.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2

   common counts, y, dy

   y = double(y)
   dy = double(dy)
;	stop
   ell_vector = penalty_function_magpd(q, 'fdrelmaxwellp')

   return, ell_vector[5]

end

;-----------------------------------------------------------------------------

function penalty_gradient_magpd_drelmaxwell, q

; Based on Paul O'Brien's C Inversion Library
; This function is called by DFPMIN.

; input:
; q = [q1, q2]: parameters that are varied to minimize the penalty function
; common block input:
; y = 5-element vector of measured counts (note: not count rates) from a single magpd telescope
; dy = estimated calibration error, fractional, i.e 0.25 = 25%
; dt = accumulation or average period (e.g., dt = 60 s)
; output: sum of penalty function over all channels
   COMPILE_OPT IDL2
   
   common counts, y, dy

   y = double(y)
   dy = double(dy)

   dldq = penalty_function_magpd_gradient(q, 'fdrelmaxwellp', 'fdrelmaxwellp_gradient')
   
   return, dldq

end


;-------------------------------------------------------------------------

pro amoeba_dfpmin_support

; specific spectral functions called by AMOEBA and DFPMIN:

; penalty_maged_relmaxwell
; penalty_gradient_maged_relmaxwell

; penalty_maged_powerlaw
; penalty_gradient_maged_powerlaw

; penalty_maged_exponential
; penalty_gradient_maged_exponential

; penalty_magpd_relmaxwell
; penalty_gradient_magpd_relmaxwell

; penalty_magpd_powerlaw
; penalty_gradient_magpd_powerlaw

; penalty_magpd_exponential
; penalty_gradient_magpd_exponential
   print, 'amoeba_dfpmin_support compiled....'

end
