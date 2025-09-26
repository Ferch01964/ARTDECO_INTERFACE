

MODULE MCONSTANTS

  ! It contains constants and paths to/for the code
  ! to be changed by the user  

  IMPLICIT NONE

  CHARACTER(LEN=5), PUBLIC, PARAMETER :: version ='v 1.0'

  !------------------------------------------------------------------------------------------------
  ! Input and Output directories

  ! ---------- TO BE CHANGED, PATH to ARTDECO -----------
  ! CHARACTER (LEN=34), PUBLIC, PARAMETER :: artdeco_path    ='/home/mcharek/Desktop/ARTDECO_REV_160/'

  !!! Line added on 10/03/14
  CHARACTER(LEN=500)  :: artdeco_path

  CHARACTER (LEN=6), PUBLIC, PARAMETER :: dir_input       ='input/'
  CHARACTER (LEN=12), PUBLIC, PARAMETER :: dir_vdist      ='input/vdist/'
  CHARACTER (LEN=9), PUBLIC, PARAMETER :: dir_prop        ='lib/prop/'
  CHARACTER (LEN=8), PUBLIC, PARAMETER :: dir_atm         ='lib/atm/'
  CHARACTER (LEN=8), PUBLIC, PARAMETER :: dir_opt         ='lib/opt/'
  CHARACTER (LEN=11), PUBLIC, PARAMETER :: dir_refind     ='lib/refind/'
  CHARACTER (LEN=12), PUBLIC, PARAMETER :: dir_surf       ='lib/surface/'
  CHARACTER (LEN=11), PUBLIC, PARAMETER :: dir_solrad     ='lib/solrad/'
  CHARACTER (LEN=9), PUBLIC, PARAMETER :: dir_kdis        ='lib/kdis/'
  CHARACTER (LEN=8), PUBLIC, PARAMETER :: dir_libihm      ='lib/ihm/'
  CHARACTER (LEN=11), PUBLIC, PARAMETER :: dir_filter     ='lib/filter/'
  CHARACTER (LEN=9), PUBLIC, PARAMETER :: dir_baum        ='lib/baum/'
  CHARACTER (*), PUBLIC, PARAMETER :: dir_betal       ='lib/betal/'


  CHARACTER (LEN=4), PUBLIC, PARAMETER :: dir_out         ='out/'

  !------------------------------------------------------------------------------------------------
  ! Fortran language related constants

  ! "kind PARAMETER"  for double precision
  INTEGER, PUBLIC, PARAMETER        :: dp        = SELECTED_REAL_KIND(P=15)
  REAL(KIND=dp), PUBLIC, PARAMETER  :: tiniestdp = tiny(1.0_dp)
  REAL(KIND=dp), PUBLIC, PARAMETER  :: hugestdp  = huge(1.0_dp)

  CHARACTER (LEN=7), PUBLIC, PARAMETER :: undef_c    = 'unknown'
  INTEGER, PUBLIC, PARAMETER       :: undef_i    = -32768
  REAL(KIND=dp), PUBLIC, PARAMETER :: undef_dp   = -32768.0_dp

  INTEGER, PUBLIC, PARAMETER :: max_len = 500

  !------------------------------------------------------------------------------------------------
  ! Mathematical and Physical constants
  REAL(KIND=dp), PUBLIC, PARAMETER :: xpi     = 3.1415926535897932384_dp ! pi (!)
  REAL(KIND=dp), PUBLIC, PARAMETER :: deg2rad = xpi / 180.0_dp


  ! line added on 13/05/15 to add Boltzmann constant
  REAL(KIND=dp), PUBLIC, PARAMETER ::  Kb  = 1.3806488D-23                ! in [m2 kg s-2 K-1]


  !REAL(KIND=dp), PUBLIC, PARAMETER :: plkavg_res = 1.0d-2

  !REAL(KIND=dp), PUBLIC :: mm_air     = 28.97_dp
  ! REAL(KIND=dp), PUBLIC :: mm_o3      = 48.00_dp
  ! REAL(KIND=dp), PUBLIC :: mm_h2o     = 18.00_dp
  ! REAL(KIND=dp), PUBLIC :: mm_co2     = 44.00_dp
  !REAL(KIND=dp), PUBLIC :: amu        = 1.660538921e-24_dp

  !REAL(KIND=dp), PUBLIC :: cp_air = 1.012_dp ! J g-1 K-1

  !!!!! Line added on 10/03/14
  REAL(KIND=dp), PUBLIC :: plkavg_res
  REAL(KIND=dp), PUBLIC :: mm_air
  REAL(KIND=dp), PUBLIC :: mm_o3
  REAL(KIND=dp), PUBLIC :: mm_h2o
  REAL(KIND=dp), PUBLIC :: mm_co2
  REAL(KIND=dp), PUBLIC :: amu
  REAL(KIND=dp), PUBLIC :: cp_air


  !------------------------------------------------------------------------------------------------
  ! code related constants

  REAL(KIND=dp), PUBLIC, PARAMETER :: EPS_LAMB  = 1e-5_dp              ! (in micron) if the difference between two 
  !                                                                        wavelengths is smaller than that, they are considered equal
  REAL(KIND=dp), PUBLIC, PARAMETER :: EPS_BETAL = 1e-3_dp              ! To test betal_0^11 (that should be 1)

  !INTEGER, PUBLIC, PARAMETER :: nang_min_betal_exp  = 1800             ! minimum number of gauss point to interpolate the 
  !                                                                      phase matrix on for the Betal expansion
  INTEGER, PUBLIC, PARAMETER :: nang_min_betal_exp  = 1800

  ! altitude grid used in computing particles distribution 
  INTEGER, PUBLIC, PARAMETER :: max_user_vdist_nalt  = 1000            ! maximum allowed number of altitudes in user particles distribution 

  INTEGER, PUBLIC, PARAMETER :: nalt_distrib  = 1000                   ! number of altitudes in the grid used for the particles distribution (log grid)
  REAL(KIND=dp), PUBLIC, PARAMETER :: min_non0_alt = 0.001_dp          ! (km) mimimum non 0 altitude used in the alt grid for the particles distribution
  !                                                                      (different from the altitude grid read in atm definition file)

  !------------------------------------------------------------------------------------------------
  ! PARTICLE OPTICAL PROPERTIES
  INTEGER, PUBLIC, PARAMETER :: ptcle_nopt = 03     ! number of optical properties per particle stored in ptcle_opt(nptcle, ptcle_nopt) 
  !---------------
  ! MIE3 constants 
  INTEGER, PUBLIC, PARAMETER :: mie3_npar   = 08    ! number of mie3 parameters to be read in prop_*.dat files
  !---------------
  ! HM constants 
  INTEGER, PUBLIC, PARAMETER :: hm_nmaxpar   = 08   ! maximum number of parameters to be read in prop_*.dat files for hm (phm, rhm or ihm)
  INTEGER, PUBLIC, PARAMETER :: phm_npar     = 02   ! number of parameters to be read in prop_*.dat files for phm
  INTEGER, PUBLIC, PARAMETER :: rhm_npar     = 03   ! number of parameters to be read in prop_*.dat files for rhm
  INTEGER, PUBLIC, PARAMETER :: ihm_npar     = 05   ! number of parameters to be read in prop_*.dat files for ihm
  REAL(kind=dp), PUBLIC, PARAMETER :: min_size_param_hm = 100.0_dp ! minimum size parameter for the ray tracing to be valid

  !------------------------------------------------------------------------------------------------
  ! Truncation
  REAL(kind=dp), PUBLIC, PARAMETER :: coeff_trunc_min = 0.0_dp ! limit for a WARNING if trunccoeff < coeff_trunc_min

  !------------------------------------------------------------------------------------------------
  ! Definition of derived types:

  TYPE, PUBLIC :: FICH    ! structure CHARACTERizing a file
     CHARACTER(LEN=max_len) :: name         ! file name
     INTEGER                :: unit         ! file unit
  END TYPE FICH

END MODULE MCONSTANTS
