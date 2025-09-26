############################################################################## 
##                                                                          ## 
## Additional configuration file used to define specific parameters         ## 
## for the surface configuraiton in ARTDECO.                                ## 
##                                                                          ## 
############################################################################## 
 
## Surface type ? 
##      - "lambert". 
##      - "brdf". 
Surface_type = "lambert" 
 
## Surface name ? 
##      - "cste" for a constant lambertian albedo. 
##      - "ocean" for an ocean brdf surface. 
##      - "<name>", the properties will then be read in the file 
##               "./lib/surface/surfalb_<name>.dat" for a lambertian surface or 
##               "./lib/surface/surfbrdf_<name>.dat" for a brdf surface. 
Surface_nature = "cste" 
 
## In the case of a lambertian surface with a constant albedo. 
Surface_albedo = 0.0173
 
## Surface temperature (K). 
## -1 to take the same temperature as the lowest layer  
##  in the atmospheric profile. 
Surface_temperature = 296.61
 
############################################# 
## Specific parameters for ocean surface.  ## 
############################################# 
 
## Wind speed in m/s. 
Wind_speed = 5.0 
 
## Salinity in ppt. 
Ocean_salinity = 34.3 
 
## Pigment concentration in mg/m^3. 
Pigment_concentration = 0.0 
 
## Is shadowing effect taken into account ? True or False. 
Shadowing = False 
 
########################################################################## 
##                                                                      ## 
## No user changes should be necessary under this line.                 ## 
##                                                                      ## 
########################################################################## 
def BasicParameters(): 
    Surface_parameters = [Surface_type, Surface_nature] 
    if Surface_nature == "cste": 
        Surface_parameters.append(float(Surface_albedo)) 
     
    return Surface_parameters, Surface_temperature 
 
def OceanParameters(): 
    return Wind_speed, Ocean_salinity, Pigment_concentration, Shadowing 
 
 
