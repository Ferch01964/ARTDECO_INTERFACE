##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for kdis mode of ARTDECO.                                                ##
##                                                                          ##
##############################################################################

## Name of the k-distribution parameterization. The program will
## then search for the file "./lib/kdis/reference/kdis_<name>.dat".
K_distribution_parameterization = "gamesw_1"

## List of gas taken into account for the gaseous absorption.
## If no gas are taken into account, Gas_list = []
Gas_list = ["co2", "o3", "h2o"]
## List of gas continuum taken into account.
## If no continuum is taken into account, Gas_continuum_list = []
Gas_continuum_list = ["co2", "o3", "h2o"]

## Instrumental filters list. The program will search for the files
## "./lib/filter/filter_<name>.dat".
## example : Instrumental_filters_list = ["parasol565", "parasol670"]
Instrumental_filters_list = []
    
## Wavelength unit might be :
##         - "cm^-1" or "cm-1".
##         - "microm", "micrometers", "mum".
## Spectral domain limits "start" and "end" may be unordered, the program
## will deal with it.
Wavelength_unit = "mum"
Wavelength_start = 0.2
Wavelength_end = 4
    
## Depolarization value :
##         - "default" : value = 0.0279 for every wavelength.
##         - user-defined value (float or int).
Depolarization_value = "default"


##########################################################################
##                                                                      ##
## No user changes should be necessary under this line.                 ##
##                                                                      ##
##########################################################################
def BasicParameters():
    global Wavelength_start
    global Wavelength_end
    global Depolarization_value

    ## Converts the wavelengths in micrometers and orders the wavelength_start
    ## and wavelength_end.
    if Wavelength_unit == "cm-1" or Wavelength_unit == "cm^-1":
        Wv1 = 10000. / Wavelength_start
        Wv2 = 10000. / Wavelength_end
    else:
        Wv1 = Wavelength_start
        Wv2 = Wavelength_end      
    Wavelength_start = min(Wv1, Wv2)
    Wavelength_end = max(Wv1, Wv2)
    
    if Depolarization_value == "default":
        Depolarization_value = -1

    gl = [s.lower() for s in Gas_list]
    gcl = [s.lower() for s in Gas_continuum_list]
        
    return K_distribution_parameterization, gl, \
        gcl, Instrumental_filters_list, \
        Wavelength_unit, Wavelength_start, Wavelength_end, \
        Depolarization_value
