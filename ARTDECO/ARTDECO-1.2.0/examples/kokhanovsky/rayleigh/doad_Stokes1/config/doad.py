##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for Vesperini Adding-doubling model                                      ##
##                                                                          ##
##############################################################################

## Should the program write the downward radiances at the surface ?
Print_downward_radiance = True
    
        
###########################################################################
##                                                                       ##
## Advanced parameters.                                                  ##
##                                                                       ##
###########################################################################

## Number of Gauss points (default: 24).
N_stream = 24

## Accuracy of adding calculation epsilon (default: 1e-5).
Adding_epsilon_accuracy = 1e-5

## Absolute maximum number of Fourier terms (default: 1000).
Maximum_fourier_terms = 1000

## Minimum number of Fourier terms in the BDRF/BPDF Fourier expansion.
## (default: 10).
Minimum_fourier_terms = 10


###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################
def BasicParameters():
    return Print_downward_radiance

def AdvancedParameters():
    return N_stream, Adding_epsilon_accuracy,Maximum_fourier_terms, \
        Minimum_fourier_terms
