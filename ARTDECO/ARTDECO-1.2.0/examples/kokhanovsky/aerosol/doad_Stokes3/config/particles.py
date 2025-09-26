##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for the definition of particles and clouds in ARTDECO.                   ##
##                                                                          ##
##############################################################################

## Definition of particles parameterization.
## Pattern:
## [particle name, use_betal flag, interp. flag, H-G flag, 
##  integrated opacity at 550 nm, vertical distribution type,
##  parameters for the vertical distribution (0, 1 or 2 parameters)]
##
## particle name : 
Particle_option = [
    #["aerosols_mono_10980-10990-cm-1", True, False, False, 0.5, "layer", 2.0],
    ["kokha_aer", False, False, False, 0.3262, "layer", 1.0],
    #["aerosols_mono_18180-18190-cm-1", True, False, False, 1.0, "layer", 2.0],
    #["stephens_st2", False, False, True, 5.0, "layer", 0.05]
                   ]

## Truncation method for the phase matrix. Possibilities :
##     - "none"
##     - "potter"
##     - "dfit"
##     - "dm"
Truncation_method = "dm"

## Number of expansion coefficients.
## If Number_of_betal = -1 or Number_of_betal = "optimized", ARTDECO will 
## automatically set it to the optimized value regarding the number of 
## quadrature points (so called number of streams) for the radiative transfer.
Number_of_betal = -1

## Should the intensity correction (TMS) of Nakajima & Tanaka (1988)
## should be applied ?
## If so, the Number_of_betal should be set to "optimized".
TMS_correction = True


###########################################################################
##                                                                       ##
## Advanced parameters for the Potter truncation method.                 ##
##  (Only used if Truncation_method = "potter")                          ##
##                                                                       ##
###########################################################################

## Theta_min (default = 14.0)
Theta_min = 14.0

## Theta_max (default = 15.0)
Theta_max = 15.0


###########################################################################
##                                                                       ##
## Advanced parameters for the d-fit truncation method.                  ##
##  (Only used if Truncation_method = "dfit")                            ##
##                                                                       ##
###########################################################################

## Theta_cut (default = 0.0)
Theta_cut = 0.0

## Fit all (default = False) ?
Fit_all = False


###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################

def AdvancedParameters():
    return Theta_min, Theta_max

##############################################################################
##                                                                          ##
## No user changes should be necessary under this line.                     ##
##                                                                          ##
##############################################################################
def BasicParameters():
    if Number_of_betal == "optimized":
        Nbl = -1
    else:
        Nbl = Number_of_betal
    
    return Particle_option, Truncation_method, Nbl, TMS_correction

def PotterParameters():
    return Theta_min, Theta_max

def DfitParameters():
    return Theta_cut, Fit_all
 
