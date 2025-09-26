

MODULE MCOMMON
 
  ! This module will contain declaration of variables to be
  ! seen all over ARTDECO

  USE MCONSTANTS, only : dp, max_len

  IMPLICIT NONE

  LOGICAL, PUBLIC :: f2py_mode = .false.
  REAL(KIND=dp), PUBLIC   :: reflamb = 0.550_dp             ! default value for reference wavelengths (in micron)     

  CHARACTER(LEN=max_len), PUBLIC :: artdeco_infile_name

  !------------- General ARTDECO keywords and options ...
  LOGICAL, PUBLIC :: verbose             ! .TRUE. --> verbose on
  LOGICAL, PUBLIC :: opt_only            ! .TRUE. --> the code will only compute the optical properties (no Betal expansion, no radiative transfer)
  LOGICAL, PUBLIC :: betal_only          ! .TRUE. --> the code will compute the betal expansion ( no radiative transfer)
  LOGICAL, PUBLIC :: do_rt               ! .TRUE. --> the code will compute radiative transfer
  LOGICAL, PUBLIC :: thermal             ! .TRUE. --> the radiative transfer will be computed accounting for in-situ thermal emission
  LOGICAL, PUBLIC :: sinsca_corint       ! .TRUE. --> single scattering intensity correction will be performed
  LOGICAL, PUBLIC :: warning             ! .TRUE. --> WARNING message are displayed
  LOGICAL, PUBLIC :: no_rayleigh 
  LOGICAL, PUBLIC :: thermal_only
  LOGICAL, PUBLIC :: flux_only

  !! line added on 12/05/15 for atmospheric profile with u_atm in ppmv
  LOGICAL, PUBLIC :: atm_in_ppmv

  LOGICAL, PUBLIC :: od_no_check  !.TRUE. --> DISORT does not perform the CHECKIN() routine

  LOGICAL, PUBLIC :: print_recomp ! for recomposed phase function to be printed out
  LOGICAL, PUBLIC :: print_betal  ! for Legendre polynomial coeff. to be printed out
  LOGICAL, PUBLIC :: print_aitaui
  LOGICAL, PUBLIC :: print_down_rad  ! for Downward intensities to be printed out

  !!!line added on 08/04/2015 for user iterface in choosing output file elements
!  LOGICAL, PUBLIC                 :: Change_Output_flux_only,Change_Output
!  INTEGER,PUBLIC                  :: alt_level(5)
!  CHARACTER(LEN=max_len), PUBLIC  :: choix
!  REAL(kind=dp), PUBLIC           :: sza_user, vza_user, vaa_user
!  INTEGER,  PUBLIC                :: isza, ivza, ivaa


  CHARACTER(LEN=max_len), PUBLIC :: mode         ! monochromatic or k-distribution
  CHARACTER(LEN=max_len), PUBLIC :: rt_model     ! the radiative transfer model to be used
  CHARACTER(LEN=max_len), PUBLIC :: trunc_method ! truncation method to be used
  CHARACTER(LEN=max_len), PUBLIC :: out_root     ! root name for output files


  !!! Line added on  20/03/14 to  directly read the name of the input file
  !!! "artdeco_in.dat"
  CHARACTER(LEN=max_len), PUBLIC :: in_root     ! root name for output files


  integer, PUBLIC :: start_date(8)
  integer, PUBLIC :: end_date(8)
  REAL, PUBLIC    :: start_tcpu
  REAL, PUBLIC    :: end_tcpu

  CHARACTER (LEN=max_len), PUBLIC :: dir_root_out

  !------------- timing variables
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: rt_cputime(:) ! rt_cputime(nlambda) CPU time while calling RT subroutine 
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: rt_cputime_filter(:) ! rt_cputime(filter) CPU time while calling RT subroutine 

  !------------- Geometrical conditions
  INTEGER, PUBLIC :: nsza                                ! number of solar zenith angle(s)
  INTEGER, PUBLIC :: nvza                                ! number of view zenith angle(s)

  INTEGER, PUBLIC :: nvza_up                             ! number of view zenith angle(s) in upward direction
  INTEGER, PUBLIC :: nvza_dw                             ! number of view zenith angle(s) in downward direction
  INTEGER, PUBLIC, ALLOCATABLE       :: ind_vza_up(:)    ! Indices for which UMU = cos(VZA) is positive
  INTEGER, PUBLIC, ALLOCATABLE       :: ind_vza_dw(:)    ! Indices for which UMU = cos(VZA) is negative

  REAL(KIND=dp), PUBLIC,  ALLOCATABLE :: temp_vza(:)
    
  INTEGER, PUBLIC :: nvaa                                ! number of view azimuthal angle(s)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: sza(:)           ! sza(nsza) : solar zenith angle(s) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: vza(:)           ! vza(nvza) : view zenith angle(s) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: vaa(:)           ! vaa(nvaa) : view azimuthal angle(s)

  ! Location and time for solar geometry computation
  LOGICAL :: get_possol ! flag to known wether we compute solar position (and constant) 
  !                       depending on time and geographical coordinate
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: lon(:)  ! lon(nsza) longitude of the scene (decimal degree)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: lat(:)  ! lat(nsza) latitude of the scene (decimal degree)
  INTEGER, PUBLIC, ALLOCATABLE       :: day(:)  ! day(nsza) day of the year 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: h_tu(:) ! h_tu(nsza) universal time
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: varsol_fact(:) ! varsol_fact(nsza) : scaling factor accounting for
  !                                                        Earth-Sun distance variation within a year.

  !------------ filters
  LOGICAL, PUBLIC :: filter
  INTEGER, PUBLIC :: nfilter
  CHARACTER (LEN=max_len), ALLOCATABLE, PUBLIC :: filter_name(:) ! filter_name(nfilter)
  INTEGER, PUBLIC :: nmaxlamb_filter ! dimensionning variable
  INTEGER, PUBLIC, ALLOCATABLE :: nlamb_filter(:)         !  nlamb_filter(nfilter)
  REAL(KIND=DP), PUBLIC, ALLOCATABLE :: bound_lamb_filter(:,:) ! bound_lamb_filter(2,nfilter) in microns
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: lamb_filter(:,:)       ! lamb_filter(nfilter, nmaxlamb_filter) in microns
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: trans_filter(:,:)      ! trans_filter(nfilter, nmaxlamb_filter)

  !------------- Incident radiation (solar beam)
  !!! line added on 12/09/14 to read Fbeam when entering by user
  LOGICAL, PUBLIC  :: Fbeam_user_flag            ! whether user enter or not the value of fbeam
  REAL(kind=dp), PUBLIC :: Fbeam_user_value      ! Value of Fbeam as constant

  CHARACTER(LEN=max_len), PUBLIC :: solrad ! monochromatic solar radiation spectrum to be used (read in file)
  INTEGER, PUBLIC :: solrad_nwvl           ! number of wavelength in solrad
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: solrad_wvl(:)  ! solrad_wvl(solrad_nwvl) solrad wavelengths (micron)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: solrad_F0(:)   ! solrad_F0(solrad_nwvl) monochromatic values solrad spectrum (W m-2 microns-1)
  REAL(kind=dp), PUBLIC  :: solrad_solcste   ! solar constant (W m-2) 

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: F0(:)  ! Solar radiation used by RTE solver F0(nlambda)
  REAL(kind=dp), PUBLIC :: Svin(4)             ! NB : the normailzed stokes vector (S_1 = 1) for 
  !                                                   incoming light is not dependent on wlambda
  !                                                   for solar light --> always {1,0,0,0}        

  !------------- Atmosphere model (Temp., Pressure, gaz concentration, etc)
  CHARACTER(LEN=max_len), PUBLIC :: atm_type                ! type for the atmosphere model (linked to a file to be read)
  INTEGER, PUBLIC :: nalt_atm                               ! number of altitude for the atmosphere definition
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: alt_atm(:)          ! alt_atm(nalt_atm)   : altitude grid (km)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: t_atm(:)            ! t_atm(nalt_atm)     : temperature (K)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: p_atm(:)            ! p_atm(nalt_atm)     : pressure (mb)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: u_air_atm(:)        ! air concentration profil (cm-3)
  INTEGER, PUBLIC   :: ngas_u_atm                              ! number of gas 
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: gas_u_atm(:)  ! (ngas_u_atm) name of gases 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: u_atm(:,:)             ! (ngas_u_atm, nalt_atm) concentration profil (cm-3)


  !--- line added on 18/09/14 to read list of gazes given in artdeco_in.dat
  INTEGER, PUBLIC :: n_list_gas                                     ! number of gas given by user
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: list_gas(:)           ! list_gas(n_list_gas) : name of gases


  ! ------------ k-distribution 
  CHARACTER(LEN=max_len), PUBLIC :: kdis_model ! the root name of kdis_ files to be read
  INTEGER, PUBLIC :: kdis_nmaxai   ! maximum number of quadrature points over the 
  !                                  different wavelengths and species
  INTEGER, PUBLIC :: kdis_nt   ! number of temperatures defined in the kdis model
  INTEGER, PUBLIC :: kdis_np   ! number of pressures defined in the kdis model
  INTEGER, PUBLIC :: kdis_nc   ! number of concentrations defined in the kdis model
  INTEGER, PUBLIC :: kdis_nwvl ! number of wavelengths defined in the kdis model
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_t(:)             ! temperatures defined in the kdis model (K)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_p(:)             ! pressures defined in the kdis model (mb)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_c(:)             ! concentrations defined in the kdis model (cm-3)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_wvlband(:,:)     ! kdis wavelengths, kdis_wvlband(3,kdis_nwvl) 
  !                                                             kdis_wvlband(1, :) central wavelength of the band
  !                                                             kdis_wvlband(2, :) lower boundary wavelength of the band
  !                                                             kdis_wvlband(3, :) higher boundary  wavelength of the band
  ! species without concentration dependency
  INTEGER, PUBLIC :: kdis_nsp  ! number of gas species accounted for in the kdis method 
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: kdis_species(:) ! gas species described in the kdis model
  INTEGER, PUBLIC, ALLOCATABLE       :: kdis_nai(:,:)         !  kdis_nai(kdis_nsp,kdis_nwvl) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_ai(:,:,:)        !  kdis_ai(kdis_nsp,kdis_nwvl,kdis_nmaxia) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_ki(:,:,:,:,:)    !  kdis_ki(kdis_nsp,kdis_nwvl,kdis_nmaxia,kdis_np,kdis_nt) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_xsect(:,:)       !  kdis_xsect(kdis_nsp,kdis_nwvl) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_fcont(:)         !  a weight to the continuum to be added 
  ! species with concentration dependency
  INTEGER, PUBLIC :: kdis_nsp_c  ! number of gas species accounted for in the kdis method
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: kdis_species_c(:) ! gas species described in the kdis model
  INTEGER, PUBLIC, ALLOCATABLE       :: kdis_nai_c(:,:)         !  kdis_nai(kdis_nsp,kdis_nwvl) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_ai_c(:,:,:)        !  kdis_ai(kdis_nsp,kdis_nwvl,kdis_nmaxia) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_ki_c(:,:,:,:,:,:)  !  kdis_ki(kdis_nsp,kdis_nwvl,kdis_nmaxia,kdis_np,kdis_nt,kdis_nc) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_xsect_c(:,:)       !  kdis_xsect(kdis_nsp,kdis_nwvl) 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: kdis_fcont_c(:)         !  a weight to the continuum to be added 

  ! solar radiation integrated over k-dis bands
  LOGICAL, PUBLIC :: kdis_solrad ! whether the integrated solar radiation file already exists for the k-dis
  INTEGER, PUBLIC :: kdis_solrad_nwvl                       ! number of wavelength 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: kdis_solrad_wvl(:)  ! wavelengths (micron)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: kdis_solrad_F0(:)   ! value of the integrated solar radiation (W m-2)
  REAL(kind=dp), PUBLIC  :: kdis_solrad_solcste             ! solar constant (W m-2)

  !------------- Working wavelength
  INTEGER, PUBLIC                    :: nlambda              ! number of wavelengths
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: wlambda(:)           ! wavelengths to work at (microns)
  INTEGER, PUBLIC, ALLOCATABLE :: lambda_ikdis(:)            ! index of kdis_wvlband(1,:) corresponding to the working wavelength
  INTEGER, PUBLIC                    :: ireflamb             ! index of the reference wavelength among working wavelength
  !                                                            ireflamb = 0 if reference wavelength does no match any working wavelength
  !                                                                         or if no RT will be performed
  REAL(KIND=dp), PUBLIC              :: lambda_bound(2)      ! limit wavelengths to work with (used in kdis mode to integrate the flux)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: depol(:)             ! depol(nlambda) : depolarization factor

  !------------- particles description

  INTEGER, PUBLIC :: nptcle                                              ! number of particle population 
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: ptcle_type(:)           ! ptcle_type(nptcle)              : particle type for each population

  LOGICAL, PUBLIC :: hg_flag                                             ! whether H-G approx. will be used for some ptcle for the RT
  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_hg(:)                            ! ptcle_hg(nptcle)          : whether we approximate the phase function with an H-G 

  !!! Line added on 19/02/14 for ptcle_use_betal
  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_use_betal(:)                            ! ptcle_use_betal(nptcle)          : whether we read an exisiting betal from file produced by user


  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_opt_interp(:)                    ! ptcle_opt_interp(nptcle)  : whether we interpolate the ptcle optical properties using stored one

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_tau_reflamb(:)             ! ptcle_tau_reflamb(nptcle)       : optical depth for each population at reference wavelength
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: ptcle_vdist_type(:)     ! ptcle_vdist_type(nptcle)        : vertical distribution type for each population 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_h_vdist(:)                 ! ptcle_h_vdist(nptcle)           : reference altitude 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_sd_gauss_vdist(:)          ! ptcle_sd_gauss_vdist(nptcle)    : standard deviation (when gaussian vertical distribution)
  INTEGER, public, allocatable :: ptcle_user_vdist_nalt(:)               ! ptcle_user_vdist_nalt(nptcle) : number of sampled altitude in user particle vertical distribution
  REAL(kind=dp), public, allocatable :: ptcle_user_vdist_alt(:,:)        ! ptcle_user_vdist_alt(nptcle,max_user_vdist_nalt) : altitude in user particle vertical distribution (km)
  REAL(kind=dp), public, allocatable :: ptcle_user_vdist(:,:)            ! ptcle_user_vdist(nptcle,max_user_vdist_nalt) : user particle vertical distribution (arbitrary units)

  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: ptcle_opt_model(:)      ! ptcle_opt_model(nptcle)         : code to be used to compute opt. prop. for each particle (if opt_exist .eqv. .false.)
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: ptcle_material(:)       ! ptcle_material(nptcle)          : particle material for each population (attached to a refractive index) 
  logical, public :: flag_mie ! whether new opt. prop. will be computed using mie
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_param_mie3(:,:)            ! ptcle_param_mie3(nptcle, mie3_npar) : mie3 parameters to be read in prop_*.dat files  
  logical, public :: flag_hm ! whether new opt. prop. will be computed using HM
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nphot_hm(:)                      !  ptcle_nphot_hm(nptcle)       : number of photons to compute HM particle optical properties
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_param_hm(:,:)              ! ptcle_param_hm(nptcle, hm_nmaxpar)  : parameters to be read in prop_*.dat files for 
  !                                                                                                              hexagonal monocrystal properties computation (phm, rhm or ihm)
  logical, public :: flag_baum ! whether new opt. prop. will be computed using Baum
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_param_baum(:)              ! ptcle_param_baum(nptcle)            : only on parameter for Baum cirrus model : Deff (microns)


  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_opt_exist_reflamb(:)             ! ptcle_opt_exist_reflamb(nptcle) : logical= .true. if opt_type_refwave.dat already exist
  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_opt_exist(:,:)                   ! ptcle_opt_exist(nptcle,nlambda) : logical= .true. if opt_type.dat already exist

  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nelem(:)                         ! ptcle_nelem(nptcle) : number of necessary element for the phase pmatrix
  !                                                                             nelem = 4 for spherical particles
  !                                                                             nelem = 6 otherwise


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!   file modified on  19/02/14      !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! BETAL_STORE : these are optical properties read in betal_[].dat files
  LOGICAL, PUBLIC :: flag_store_betal
  INTEGER, PUBLIC :: nbetal_max_store
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_opt_betal_store(:,:,:)           ! ptcle_betal_store(nptcle, ptcle_nopt, nlamb_max_store)  : betal for particles optical properties 
  !                                                                          ptcle_betal_store(nptcle, 1,lamb) : Cext   : the extinction coeff 
  !                                                                          ptcle_betal_store(nptcle, 2,lamb) : w0       : the single scattering albedo
  !                                                                          ptcle_betal_store(nptcle, 3,lamb) : g      : asymetry factor
 INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nlamb_betal_store(:)             ! ptcle_nlamb_betal_store(nptcle) : number of wavelengths

 INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nbetal_betal_store(:,:)          ! ptcle_nbetal_betal_store(nptcle,nlamb_max_store) : number of betal 
 REAL, PUBLIC, ALLOCATABLE ::ptcle_trunc_coeff_betal_store(:,:)          ! ptcle_trunc_coeff_betal_store(nptcle,nlamb_max_store) : truncation coefficient

 REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_lamb_betal_store(:,:)      ! ptcle_lamb_phasemat_store(nptcle, nlamb_max_store)  wavelengths 
 REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_betal_store(:,:,:,:)       ! ptcle_betal_store(nptcle, 4, nang_max_store, nlamb_max_store) : the betal coeffs
  !                                                                          alpha1 =  ptcle_betal_store(:, 1, :, :)
  !                                                                          alpha2 =  ptcle_betal_store(:, 2, :, :)
  !                                                                          alpha3 =  ptcle_betal_store(:, 3, :, :)
  !                                                                          beta1 =   ptcle_betal_store(:, 4, :, :)




  ! OPT_STORE : these are optical properties read in opt_[].dat files
  LOGICAL, PUBLIC :: flag_store_opt
  INTEGER, PUBLIC :: nlamb_max_store
  INTEGER, PUBLIC :: nang_max_store 
  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_flag_store_opt(:)                ! ptcle_flag_store_opt(nptcle) :  whether there are stored optical properties for each ptcle
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_opt_store(:,:,:)           ! ptcle_opt_store(nptcle, ptcle_nopt, nlamb_max_store)  : particles optical properties 
  !                                                                              ptcle_opt_store(nptcle, 1,lamb) : Cext   : the average extinction cross section 
  !                                                                              ptcle_opt_store(nptcle, 2,lamb) : w0     : the single scattering albedo
  !                                                                              ptcle_opt_store(nptcle, 3,lamb) : g      : asymetry factor
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nlamb_phasemat_store(:)          ! ptcle_nlamb_phasemat_store(nptcle) : number of wavelengths
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nang_phasemat_store(:,:)         ! ptcle_nang_phasemat_store(nptcle,nlamb_max_store) : number of angle the phase matrix is defined on
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_lamb_phasemat_store(:,:)   ! ptcle_lamb_phasemat_store(nptcle, nlamb_max_store)  wavelengths
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_u_phasemat_store(:,:,:)    ! ptcle_u_phasemat_store(nptcle, nang_max_store, nlamb_max_store)  : cos(theta) 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat_store(:,:,:,:)    ! ptcle_phasemat_store(nptcle, 6, nang_max_store, nlamb_max_store) : the phase matrix
  !                                                                          F11 =  ptcle_phasemat_store(:, 1, :, :)
  !                                                                          F22 =  ptcle_phasemat_store(:, 2, :, :)
  !                                                                          F33 =  ptcle_phasemat_store(:, 3, :, :)
  !                                                                          F44 =  ptcle_phasemat_store(:, 4, :, :)
  !                                                                          F21 =  ptcle_phasemat_store(:, 5, :, :)
  !                                                                          F34 =  ptcle_phasemat_store(:, 6, :, :)

  ! OPT_NEW : these are optical properties computed during the run
  LOGICAL, PUBLIC :: flag_new_opt
  INTEGER, PUBLIC :: nlamb_max_new
  INTEGER, PUBLIC :: nang_max_new 
  LOGICAL, PUBLIC, ALLOCATABLE :: ptcle_flag_new_opt(:)                ! ptcle_flag_new_opt(nptcle) :  whether there are new optical properties for each ptcle
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_opt_new(:,:,:)           ! ptcle_opt_new(nptcle, ptcle_nopt, nlamb_max_new)  : particles optical properties 
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nlamb_phasemat_new(:)          ! ptcle_nlamb_phasemat_new(nptcle) : number of wavelengths
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nang_phasemat_new(:,:)         ! ptcle_nang_phasemat_new(nptcle,nlamb_max_new) : number of angle the phase matrix is defined on
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_lamb_phasemat_new(:,:)   ! ptcle_lamb_phasemat_new(nptcle, nlamb_max_new)  wavelengths
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_u_phasemat_new(:,:,:)    ! ptcle_u_phasemat_new(nptcle, nang_max_new, nlamb_max_new)  : cos(theta) 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat_new(:,:,:,:)    ! ptcle_phasemat_new(nptcle, 6, nang_max_new, nlamb_max_new) : the phase matrix

  ! these are optical properties used to get Betal and perform RT.
  ! these are obtained from OPT_NEW or OPT_STORE
  INTEGER, PUBLIC :: nang_max
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_ilambda(:,:) ! ptcle_ilambda(nptcle, nlambda) : This variable keeps track of where to find wlambda()
  !                                                                                     in ptcle_store or ptcle_new arrays: 
  !                                                         EX : for wlambda(j) and ptcle_type(i) 
  !                                                         if ptcle_opt_exist(i,j) is TRUE you have  : ptcle_lamb_phasemat_store(i,ptcle_ilambda(i,j)) = wlambda(j)
  !                                                         if ptcle_opt_exist(i,j) is FALSE you have : ptcle_lamb_phasemat_new(i,ptcle_ilambda(i,j))   = wlambda(j)
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_ireflamb(:)  ! ptcle_ireflamb(nptcle) : This variable keeps track of where to find the reflamb
  !                                                                                     in ptcle_store or ptcle_new arrays: 
  !                                                         EX : for  ptcle_type(i) 
  !                                                         if ptcle_opt_exist_reflamb(i) is TRUE you have  : ptcle_lamb_phasemat_store(i,ptcle_ireflamb(i)) = reflamb
  !                                                         if ptcle_opt_exist_reflamb(i) is FALSE you have : ptcle_lamb_phasemat_new(i,ptcle_ireflamb(i))   = reflamb
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_Cext_reflamb(:)            ! ptcle_Cext_reflamb(nptcle) 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_opt(:,:,:)                 ! ptcle_opt(nptcle, ptcle_nopt, nlambda)  : particles optical properties 
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nang_phasemat(:,:)               ! ptcle_nang_phasemat(nptcle,nlambda) : number of angle the phase matrix is defined on
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_u_phasemat(:,:,:)          ! ptcle_u_phasemat(nptcle, nang_max, nlambda)  : cos(theta) 
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat(:,:,:,:)          ! ptcle_phasemat(nptcle, 6, nang_max, nlambda) : the phase matrix


  INTEGER, PUBLIC :: nbetal_in                                           ! read on artdeco_in.dat
  INTEGER, PUBLIC :: nmaxbetal                                           ! dimensioning variable : highest encountered order of Legendre expension coefficient
  INTEGER, PUBLIC, ALLOCATABLE :: ptcle_nbetal(:,:)                      ! ptcle_nbetal(nptcle,nlambda)     : highest order of Legendre expension coefficient used to represent any 
  !                                                                                                           element of the single scattering phase matrix
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_betal(:,:,:,:)             ! ptcle_betal(nptcle, 6, 0:nmaxbetal, nlambda) : Legendre coefficients for alpha1, alpha2 
  !                                                                                                              alpha3, alpha4, beta1, beta2 elements of the 
  !                                                                                                              single scattering phase matrix for each population
  !                                                                                                              alpha1 = ptcle_betal(:, 1,:,:)
  !                                                                                                              alpha2 = ptcle_betal(:, 2,:,:)
  !                                                                                                              alpha3 = ptcle_betal(:, 3,:,:)
  !                                                                                                              alpha4 = ptcle_betal(:, 4,:,:)
  !                                                                                                              beta1  = ptcle_betal(:, 5,:,:)
  !                                                                                                              beta2  = ptcle_betal(:, 6,:,:)

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_trunc_normphasemat(:,:)    ! ptcle_trunc_normphasemat(nptcle, nlambda)  : the integral of the recomposed phase function (should be 2)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_trunccoeff(:,:)            ! ptcle_trunccoeff(nptcle,nlambda) : truncation coefficient
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_trunc_phasemat(:,:,:,:)    ! ptcle_trunc_phasemat(nptcle, 6, nang_max,nlambda) : phase matrix recomposed after truncation

  ! material refractive index read in files
  INTEGER, PUBLIC :: nmaterial                                              ! number of different materials that is needed to compute optical properties
  CHARACTER(LEN=max_len), PUBLIC, ALLOCATABLE :: material_name(:)           ! material(nmaterial)              : material name corresponding to material_refind
  INTEGER, PUBLIC  :: nlambmax_refind                                       ! dimensioning varible:  maximum allowed number of wavelenth for the 
  !                                                                           refractive index arrays to be read in files
  INTEGER, PUBLIC, ALLOCATABLE :: material_nlamb_refind(:)                  ! material_nlamb_refind(nmaterial) : number of wavelength
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: material_refind(:,:,:)              ! material_refind(nmaterial,3,nlambmax_refind) refractive index read in files for used material
  !                                                                               material_refind(:,1, :) : wavelength in microns
  !                                                                               material_refind(:,2, :) : real part
  !                                                                               material_refind(:,3, :) : imaginary part 

  !------------ Layering of the atmosphere (for RTE solvers) -----------------------------------

  INTEGER, PUBLIC :: nlayers        ! number of layers to discretize the atmosphere into 
  !                                   for radiative transfer computation (call of rt_model)

  INTEGER, PUBLIC :: nesp_layers   ! = kdis_nsp + kdis_nspc_c : total number of absorbing species
  !                                   NOTE : if mode .ne. 'kdis', nesp_layers = 1
  INTEGER, PUBLIC :: naimax_layers ! = kdis_nmaxai, maximum number of quadrature points over the 
  !                                                 different wavelengths and species
  !                                                 NOTE1 : it is a dimensioning variable
  !                                                 NOTE2 : if mode .ne. 'kdis' naimax_layers = 1
  INTEGER, PUBLIC, ALLOCATABLE :: nai_layers(:,:)  ! nai_layers(nesp_layers, nlambda) : number of quadrature points for all the 
  !                                                                                     different wavelengths and species
  !                                                                                     if mode .ne. 'kdis' nai_layers(:,:) = 1
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ai_layers(:,:,:)              ! ai_layers(nesp_layers, naimax_layers, nlambda) : kdis weight 


  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: gas_dtauabs_layers(:,:,:,:)    ! gas_dtauabs_layers(nlayers, nesp_layers, naimax_layers, nlambda)   : 
  !                                                                                   gas absorption optical thickness for each layer, each k-dis quadrature 
  !                                                                                   point of each absorbing species and at each wavelength

  ! For DISORT and ADDING codes
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_dtauabs_layers(:,:) ! ptcle_dtauabs_layers(nlayers,nlambda) : absorption optical thickness due to particles

  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: tot_dtausca_layers(:,:)   ! tot_dtausca_layers(nlayers, nlambda) : Total (gas+particles) optical scattering optical thickness 
  INTEGER, PUBLIC, ALLOCATABLE :: nbetal_layers(:)                ! nbetal_layers(nlambda) : highest order of Legendre expension coefficient used to represent any 
  !                                                                                          element of the single scattering phase matrix for a given layer
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: betal_layers(:,:,:,:)     ! betal_layers(nlayers,6,0:nmaxbetal,nlambda) : Legendre coefficients for alpha1, alpha2 
  !                                                                                                           alpha3, alpha4, beta1, beta2 elements of the 
  !                                                                                                           single scattering phase matrix for each layers 

  ! in case of MC or single scattering radiative transfer:
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_dtau_layers(:,:)         ! ptcle_dtau_layers(nlayers,nlambda) : extinction optical thickness due to particles
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_ssa_layers(:,:)          ! ptcle_ssa_layers(nlayers,nlambda)  : single scattering albedo of particles only in a layers
  INTEGER, PUBLIC, ALLOCATABLE       :: ptcle_nang_phasemat_layers(:)  ! ptcle_nang_phasemat_layers(nlambda)                          : number of angles the phase matrix is defined on
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_u_phasemat_layers(:,:)   ! ptcle_u_phasemat_layers(nang_max, nlambda) : cos(theta)
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat_layers(:,:,:,:) ! ptcle_phasemat_layers(nlayers, 6, nang_max,nlambda) : phase matrix for each layers for particles only
  !                                                                                                                                     NO RAYLEIGH CONTRIBUTION 
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: gas_dtausca_layers(:,:)        ! gas_dtausca_layers(nlayers, nlambda)    : scattering optical thickness due to gas


  ! for TMS corection
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_dtau_layers_unscaled(:,:)     ! ptcle_dtau_layers_unscaled(nlayers,nlambda)   
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_ssa_layers_unscaled(:,:)      ! ptcle_ssa_layers_unscaled(nlayers,nlambda)   
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat_tms_layers(:,:,:,:)  ! ptcle_phasemat_tms_layers(nlayers, 6, nang_max,nlambda)
  !                                                                           phase function used in TMS correction; actual phase
  !                                                                           function divided by (1-FLYR*SSALB)

  ! The following variables are monochromatic and computed when calling RTE solver
  ! It could not be set for all wavelengths because nlambda x nmax_aitaui_mono_layers can be very big
  ! To understand what "effective number of (ai, taui) resulting from the coupling of the (ai,ki) of each species" means, see eq.31 in Lacis and Oinas, 1991
  INTEGER, PUBLIC :: nmax_aitaui_mono_layers  ! is the maximum (over wavelengths) effective number of (ai, taui)
  !                                             resulting from the coupling of the (ai,ki) of each species.
  !                                             NOTE1 : it is a dimensioning variable
  !                                             NOTE2 : if mode .ne. 'kdis' nmax_aitaui=1
  INTEGER, PUBLIC :: n_aitaui_mono_layers     !  effective number of (ai, taui)
  !                                                         resulting from the coupling of the (ai,ki) of each species
  !                                                         NOTE : if mode .ne. 'kdis' n_aitaui_layers = 1
  INTEGER, PUBLIC  :: n_aitaui_total         ! the sum of n_aitaui_mono_layers over wavelengths
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ai_mono_layers(:)           ! ai_mono_layers(nmax_aitaui_mono_layers) : kdis weight 
  ! For DISORT and ADDING codes
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: tot_dtau_mono_layers(:,:)   ! tot_dtau_mono_layers(nlayers, nmax_aitaui_mono_layers) : Total optical thickness for each layer 
  !                                                                                                                            for each effective (ai, taui)
  !                                                                                                                            resulting from the coupling of the (ai,ki) of each species
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: ssa_mono_layers(:,:)         ! ssa_mono_layers(nlayers, nmax_aitaui_mono_layers) : single scattering albedo for each layers
  ! in case of MC or single scattering radiative transfer :
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: gas_dtau_mono_layers(:,:)  ! gas_dtau_mono_layers(nlayers,nmax_aitaui_mono_layers)   total gas absorption for each effective (ai, taui)
  !                                                                                                                            resulting from the coupling of the (ai,ki) of each species
  REAL(KIND=dp), PUBLIC, ALLOCATABLE :: gas_ssa_mono_layers(:,:)   ! gas_ssa_mono_layers(nlayers, nmax_aitaui_mono_layers)  gas SSA for each effective (ai, taui)
  !                                                                                                                            resulting from the coupling of the (ai,ki) of each species

  !------------ Surface

  CHARACTER(LEN=max_len), PUBLIC :: surface_family !  lambert or brdf
  CHARACTER(LEN=max_len), PUBLIC :: surface_type   !  ocean, bog, vegetation etc... 

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: surface_albedo(:) ! surface_albedo(nlambda) 
  CHARACTER(LEN=max_len), PUBLIC :: surface_brdf_model
  INTEGER, PUBLIC :: surface_n_par_brdf
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: surface_par_brdf(:,:) ! surface_par_brdf(surface_n_par_brdf, nlambda) 
  INTEGER, PUBLIC :: surface_n_par_bpdf
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: surface_par_bpdf(:,:) ! surface_par_bpdf(surface_n_par_bpdf, nlambda) 

  REAL(kind=dp), PUBLIC :: surface_temp     ! temperature of the surface in K

  ! ocean definition :
  REAL(kind=dp), PUBLIC :: wind_spd         ! wind speed in m/s for ocean surface BRDF computation
  INTEGER, PUBLIC       :: ocean_distr      ! wave facet statistical distribution
  !                                           1 : Cox & Munk anisotropic with Gram Charlier series correction term
  !                                           2 : Cox & Munk anisotropic
  !                                           3 : Cox & Munk isotropic
  LOGICAL, PUBLIC :: ocean_shadow           ! whether we apply or not shadowing due to 
  !                                           wave sun exposition
  REAL(kind=dp), PUBLIC :: ocean_xsal       ! ocean salinity (ppt)
  REAL(kind=dp), PUBLIC :: ocean_pcl        ! pigment concentration (in mg/m^3)
  REAL(kind=dp), PUBLIC :: ocean_Rsw = -1.0  ! Reflectance emerging from sea water (above surface) 
  !                                           If  = -1, computed in the sub. brdf_ocean as a function 
  !                                           of wavelength and ocean_pcl, ocean_xsal
  !                                           NB : should be wavelength dependant


  !------------ adding code specification variables
  REAL(kind=dp), PUBLIC :: ad_epsilon        ! accuracy of adding calculation epsilon
  INTEGER, PUBLIC :: ad_nmug                 ! number of Gauss points
  INTEGER, PUBLIC :: ad_mstop                ! Fourier stop for special geometries (1) or all angles (0)
  INTEGER, PUBLIC :: ad_nfoumax              ! absolute maximum number of Fourier terms

 
  !------------ DISORT code specification variables
  INTEGER, PUBLIC :: od_nstr                  ! Number of computational polar angles
  REAL(kind=dp), PUBLIC :: od_accur           ! Convergence criterion for azimuthal (Fourier cosines) series

  !!! Line added on 25/03/14
!  LOGICAL, PUBLIC  :: od_flag_user_Fbeam            ! whether user enter or not the value of fbeam
 ! REAL(kind=dp), PUBLIC :: od_Fbeam_user_value      ! Value of Fbeam as constant

  !------------ doad adding code specification variables
  REAL(kind=dp), PUBLIC :: doad_eps         ! accuracy of doad calculation epsilon
  INTEGER, PUBLIC :: doad_nmug              ! number of Gauss points
  INTEGER, PUBLIC :: doad_nfoumx            ! absolute maximum number of Fourier terms
  INTEGER, PUBLIC :: doad_ncoef_min_bdrf    ! minimum number of Fourier terms for the 
  !                                           BDRF/BPDF Fourier expansion

  !!! Line added on 25/03/14
 ! LOGICAL, PUBLIC  :: doad_flag_user_Fbeam              ! whether user enter or not the value of fbeam
 ! REAL(kind=dp), PUBLIC :: doad_Fbeam_user_value        ! Value of Fbeam as constant

  !----------- MCRAD1D specification
  INTEGER(kind=dp), PUBLIC :: mc_nphot ! number of photons to use for the Monte-Carlo simulation
  INTEGER, PUBLIC :: mc_nmu            ! number of angles to use for phase matrix
  INTEGER, PUBLIC :: mc_vrm                  ! to know whether we use VRM and which one (see Buras and Mayer, 2012, JQSRT, 112, 434)
  REAL(kind=dp), PUBLIC :: mc_epsilon_ddis   ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  INTEGER, PUBLIC :: mc_vrm_n_firstcp        ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  INTEGER, PUBLIC :: mc_vrm_n_sccp           ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  INTEGER, PUBLIC :: mc_vrm_n_lecp           ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  INTEGER, PUBLIC :: mc_nmax_interact        ! to limit the number of scattering (including surface reflection)
  REAL(kind=dp), PUBLIC :: mc_wspl_crit      ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  INTEGER, PUBLIC :: mc_nsplmax              ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  REAL(kind=dp), PUBLIC :: mc_wrr_crit       ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434
  REAL(kind=dp), PUBLIC :: mc_wrr_min        ! VRM : see Buras and Mayer, 2012, JQSRT, 112, 434

  !!! Line added on 25/03/14
!  LOGICAL, PUBLIC  :: mc_flag_user_Fbeam              ! whether user enter or not the value of fbeam
!  REAL(kind=dp), PUBLIC :: mc_Fbeam_user_value        ! Value of Fbeam as constant

  !----------- Potter truncation specification
  REAL(kind=dp), PUBLIC :: potter_theta_min ! For Potter truncation
  REAL(kind=dp), PUBLIC :: potter_theta_max ! 

  !----------- dfit truncation specification
  REAL(kind=dp), PUBLIC :: dfit_thetac ! angle below which we zero the weight for fitting
  LOGICAL, PUBLIC       :: dfit_fitall ! whether we fit Fij=F11 or just develop it

  !------------- 
  INTEGER, PUBLIC :: nmat                               ! number of elements of the Stokes vector taken into     
  !                                                       account (4 = full polarization, 3 = 3x3 approximation
  !                                                       2 = illegal, 1 = scalar)  

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR(:,:,:,:,:)     ! SvR(nsza, nvza, nvaa, nnmat, nlambda) : stokes vector for reflected intensities

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_DW(:,:,:,:,:)     ! SvR_DW(nsza, nvza, nvaa, nnmat, nlambda) : stokes vector for straight Downward reflected intensities


  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_sig(:,:,:,:,:) ! SvR_sig(nsza, nvza, nvaa, nnmat, nlambda) : stokes vector standart deviation for reflected intensities
  !                                                                                                         when RT is done with MCRAD1D
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_notms(:,:,:,:,:)       ! SvR_notms(nsza, nvza, nvaa, nnmat, nlambda)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_sig_notms(:,:,:,:,:)   ! SvR_sig_notms(nsza, nvza, nvaa, nnmat, nlambda) 

  REAL(kind=dp), PUBLIC, ALLOCATABLE :: flux_out(:,:,:,:)             ! flux_out(nsza,3,nalt_atm,nlambda) down/up/directional down flux for given altitude
  !                                                                     and for each wavelength band (W m-2)  
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: fluxint_out(:,:,:)            ! fluxint_out(nsza,3,nalt_atm)      down/up/directional down flux for given altitude 
  !                                                                     integrated over the wavelength range (W m-2) 
  !                                                                     only computed if no filters

  ! radiative values integrated over a filter transmission
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_filter(:,:,:,:,:)         ! SvR_filter(nsza, nvza, nvaa, nmat, filter)
  !                                                                     if TMS is performed, SvR_filter is obtained from
  !                                                                     the corrected radiances ( i.e. SvR(:,:,:,:,:) )
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: SvR_sig_filter(:,:,:,:,:)     ! SvR_sig_filter(nsza, nvza, nvaa, nmat, nfilter)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: flux_filter(:,:,:,:)          ! flux_filter(nsza, 3, nalt_atm, nfilter)


  !----------------------------------------
  !  Baum cirrus model to be read in sub. READ_DATA and interpolated in GET_OPT
  INTEGER, PUBLIC :: nlamb_baum
  INTEGER, PUBLIC :: nde_baum
  INTEGER, PUBLIC :: nmu_baum
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: lamb_baum(:)   ! lamb_baum(nlamb_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: de_baum(:)     ! de_baum(nde_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: mu_baum(:)     ! mu_baum(nmu_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: Cext_baum(:,:) ! Cext_baum(nlamb_baum, nde_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ssa_baum(:,:)  ! ssa_baum(nlamb_baum, nde_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: g_baum(:,:)    ! g_baum(nlamb_baum, nde_baum)
  REAL(kind=dp), PUBLIC, ALLOCATABLE :: ptcle_phasemat_baum(:,:,:,:)    ! ptcle_phasemat_baum(nlamb_baum, nde_baum, 6, nmu_baum) : the phase matrix
  !                                                                          F11 =  ptcle_phasemat_baum(:, :, 1, :)
  !                                                                          F22 =  ptcle_phasemat_baum(:, :, 2, :)
  !                                                                          F33 =  ptcle_phasemat_baum(:, :, 3, :)
  !                                                                          F44 =  ptcle_phasemat_baum(:, :, 4, :)
  !                                                                          F21 =  ptcle_phasemat_baum(:, :, 5, :)
  !                                                                          F34 =  ptcle_phasemat_baum(:, :, 6, :)


END MODULE MCOMMON
