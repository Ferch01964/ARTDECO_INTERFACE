##############################################################################
##                                                                          ##
## Additional configuration file used to define specific parameters         ##
## for 1D Monte Carlo radiation model                                       ##
##                                                                          ##
##############################################################################

        
###########################################################################
##                                                                       ##
## Advanced parameters.                                                  ##
##                                                                       ##
###########################################################################

## Number of photons (default: 1e6).
N_photons = 1.e5

## Maximum number of interactions (default: -1).
Maximum_interactions = -1

## Use of variance reduction methods (default: False).
Variance_reduction_methods = False

## Epsilon_ddis (default: 0.1).
Epsilon_ddis = 0.1

## Cloned photon specifications for N-tuple LE.
##     - mc-vrm-n-firstcp (default: 1).
##     - mc-vrm-n-sccp (default: 12).
##     - mc-vrm-n-lecp (default: 10).
Mc_vrm_n_firstcp = 1
Mc_vrm_n_sccp = 12
Mc_vrm_n_lecp = 10

## Critical weight for splitting (default: 3.0).
Critical_weight_splitting = 3.

## Maximum number of secondary photons / splitting (default: 5000)
Maximum_number_secondary_photons = 5000

## Critical weight for russian roulette (default: 0.3).
Critical_weight_russian_roulette = 0.3

## Minimum probability to kill a photon (default: 0.2).
Minimum_probability_photon = 0.2


###########################################################################
##                                                                       ##
## No user changes should be necessary under this line.                  ##
##                                                                       ##
###########################################################################
def AdvancedParameters():
    return N_photons, Maximum_interactions, Variance_reduction_methods,\
        Epsilon_ddis, Mc_vrm_n_firstcp, Mc_vrm_n_sccp, Mc_vrm_n_lecp,\
        Critical_weight_splitting, Maximum_number_secondary_photons,\
        Critical_weight_russian_roulette, Minimum_probability_photon
