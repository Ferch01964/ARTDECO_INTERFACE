###############################################################################
##                                                                           ##
## Additional configuration file used to define specific parameters          ##
## for the monochromatic mode of ARTDECO.                                    ##
##                                                                           ##
###############################################################################

## Wavelength unit might be :
##         - "cm^-1" or "cm-1".
##         - "microm", "micrometers", "mum".
## Wavelengths in the list must be in ascendant or descendant order, 
##the program will deal with it.
Wavelength_unit = "mum"
Wavelength_list = [0.412]#, 0.91]

## Depolarization mode :
##         - "default" : value = 0.0279 for every wavelength.
##         - user-defined value (float or integer)
Depolarization_value = 0.0


###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################
def BasicParameters():
    global Depolarization_value

    ## Converts the wavelengths in micrometers.
    if Wavelength_unit == "cm-1" or Wavelength_unit == "cm^-1":
        for i in range(len(Wavelength_list)):
            Wavelength_list[i] = 10000. / Wavelength_list[i]
    
    if Depolarization_value == "default":
        Depolarization_value = -1
        
    return Wavelength_unit, Wavelength_list, Depolarization_value
