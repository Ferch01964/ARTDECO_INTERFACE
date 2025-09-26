## In which mode is the code running ?
##      - "mono" for monochromatic mode.
##      - "kdis" for k-distribution parameterization.
Mode = "mono"

## Are there particles (or cloud) ?
With_particles = True

## Is there an atmosphere ?
With_atmosphere = False

## Which radiative transfer model will be used ?
##      - "none" : only optical properties will be computed (no RT model).
##      - "disort" : Discrete ordinate model DISORT 2.0.
##      - "doad" : Verperini Doubling-adding model.
##      - "mcrad1d" : MCRAD1D Monte Carlo model.
##      - "sinsca" : Single scattering approximation.
RTE_solver = "doad"

## If no RT model is used, should ARTDECO compute at least Legendre 
## polynomial expansion and phase matrix truncation ? If no, ARTDECO
## will only compute optical properties.
Legendre_and_truncation = True

## Is the Rayleigh molecular scattering taken into account ?
Rayleigh_scattering = True

## Which level of polarization should be taken into account ?
##      - 1 : scalar.
##      - 3 : I, Q, U components of Stokes vector.
##      - 4 : All components of Stokes vector.
Stokes_components = 3

#############################
## Program  options.

## Prefix used for the configuration files which will be written.
## These configuration files will be written in the directory "./input".
Prefix_of_input_files = "kokha_ref"

Verbose = True
Warnings = True
Log_file = "log_Test.txt"

#############################
## Output options.

## Name of the output subdirectory. It is located in "./out/".
## If it does not exist, it is created. 
## If it exists, files in it will be erased.
Output_directory = "kokhanovsky_aerosol_doad_Stokes3"

## Format of output files :
##      - "ascii" for ascii files only. 
##      - "netcdf" or "NetCDF" for NetCDF (not implemented yet).
##      - "hdf" or "HDF" for HDF (not implemented yet).
##      - "all" for all HDF and NetCDF (not implemented yet).
## In all cases, ascii files will be written.
Output_format = "ascii"

## Should the Legendre expansion coefficients be written in the output files ?
Print_betal = True

## Should the recomposed phase matrix from the Legendre
## coefficients be written in the output files ?
Print_recomposed_phase_matrix = True


##############################################################################
##                                                                          ##
## No user changes should be necessary under this line.                     ##
##                                                                          ##
##############################################################################
def BasicParameters():
    return Mode, With_particles, With_atmosphere, RTE_solver, \
        Legendre_and_truncation, Rayleigh_scattering, Stokes_components, \
        Prefix_of_input_files, Verbose, Warnings, Log_file, Output_directory, \
        Output_format, Print_betal, Print_recomposed_phase_matrix
