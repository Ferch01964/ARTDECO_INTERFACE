##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for the atmospheric profile.                                             ##
##                                                                          ##
##############################################################################

## Name of the atmospheric profile. The program will then be read in
## "./lib/atm/reference/atm_<name>.dat".
Atmospheric_profile = "day_ppmv"

## Are the values of sampled concentration in ppmv ?
Atmosphere_in_ppmv = True

## List of species with a uniform concentration distribution
## which should not be read in the atmospheric profile file.
## Format : [[name1, value1], [name2, value2]] or [].
Uniform_gas_list = [["co2", 360],
                    ["o2", 2.09e5]]

###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################
    
def BasicParameters(): 
    ugl = []
    for item in Uniform_gas_list:
        ugl.append([item[0].lower(), item[1]])

    ## Returns the parameters.
    return Atmosphere_in_ppmv, Atmospheric_profile, ugl

