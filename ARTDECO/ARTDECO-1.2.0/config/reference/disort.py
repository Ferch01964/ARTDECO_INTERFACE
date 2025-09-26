##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for DISORT model.                                                        ##
##                                                                          ##
##############################################################################

## Should radiances be computed ? Otherwise, only fluxes will be.
Compute_radiance = True

## Should the program write the downward radiances at the surface ?
Print_downward_radiance = False

## Is the thermal local emission taken into account ?
Thermal_emission = False
    
## Thermal local emission only ? In this case, the incoming radiation at 
## the top of atmosphere is set to zero.
Thermal_only = False
    
###########################################################################
##                                                                       ##
## Advanced parameters.                                                  ##
##                                                                       ##
###########################################################################

## Number of streams (default: 8).
N_stream = 8

## Convergence criterion (default: 0.0).
## (it should be between 0.0 and 0.001).
Convergence_criterion = 0.000


###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################

def BasicParameters():
    return Compute_radiance, Print_downward_radiance, Thermal_only, \
        Thermal_emission

def AdvancedParameters():        
    return N_stream, Convergence_criterion
