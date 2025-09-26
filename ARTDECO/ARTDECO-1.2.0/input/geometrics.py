##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for the geometric configurations in ARTDECO.                             ##
##                                                                          ##
##############################################################################

################################
## Incident radiation.        ##
################################

## In kdis mode : name of the incident spectrum characterization.
## The program will then read the file "./lib/solrad/solrad_<name>.dat"
Incident_spectrum = "kurudz_medium"

## In monochromatic mode, is the beam source specified by the user ? 
## Otherwise its value will be set to 1.
Fbeam_user = False
## If True, which value ?
Fbeam_value = 3.141592653

################################
## Solar configuration.       ##
################################

## Solar configuration mode :
##      - "angle" : Solar zenith angles will be defined.
##      - "position" : lat-lon positions and time will be defined.
Solar_configuration_mode = "angle"

## if Solar_configuration_mode = "angle", list of solar 
## zenith angles (in degree).
Solar_zenith_angle_list = [29.54]

## if Solar_configuration_mode = "position", list of position-time
## with the format :
## Geographic_position = [
##             [longitude1, latitude1, day1 of the year, time1 (decimal, UT)],
##             [longitude2, latitude2, day2 of the year, time2 (decimal, UT)],
##             ...,
##             ]
Geographic_position = [
	]

################################
## Observation geometry.      ##
################################

## List of view zenith angles (in degree).
View_zenith_angle_list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 
                          10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 
                          20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 
                          30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 
                          40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 
                          50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 
                          60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 
                          70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 
                          80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0,
			  90.0]

## List of view azimuthal angles (in degree).
View_azimuthal_angle_list = [0.0, 90.0, 180.0,]
    

    ##########################################################################
    ##                                                                      ##
    ## No user changes should be necessary under this line.                 ##
    ##                                                                      ##
    ##########################################################################
def BasicParameters(mode):
    ## FBeam specification.
    FBeam = []
    if Fbeam_user:
        FBeam.append(Fbeam_value)
    
    Solar_configuration_list = []
    if Solar_configuration_mode == "angle":
        for i in range(len(Solar_zenith_angle_list)):
            Solar_configuration_list.append([Solar_zenith_angle_list[i]])
    elif Solar_configuration_mode == "position":
        for i in range(len(Geographic_position)):
                Solar_configuration_list.append(Geographic_position[i])
    
    return Incident_spectrum, FBeam, \
        Solar_configuration_mode, Solar_configuration_list, \
        View_zenith_angle_list, View_azimuthal_angle_list
                         

    

