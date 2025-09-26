
MODULE MGET_OPT

  IMPLICIT NONE
  !suppress those attribute for f2py 
  !PUBLIC :: GET_OPT
  !PRIVATE :: GET_REFIND, CALL_MIE3, CALL_PHM, CALL_RHM, CALL_IHM, HG, GET_BAUM

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE GET_OPT

    USE MCONSTANTS, only : dp, &
         hm_nmaxpar, &
         mie3_npar, &
         ptcle_nopt, &
         max_len,&
         undef_c, &
         undef_i, &
         undef_dp
    !---- line added on 19/02/14 for ptcle_use_betal, ...
    USE MCOMMON, only : opt_only, &
         nmat, &
         betal_only, &
         verbose, &
         nlambda, &
         wlambda, &
         reflamb, &
         ireflamb, &
         nptcle, &
         hg_flag, &
         ptcle_opt_exist, &
         ptcle_opt_exist_reflamb, &
         ptcle_opt_model, &
         ptcle_type, &
         ptcle_hg, &
         ptcle_nelem, &
         ptcle_opt_interp, &
         
         ptcle_use_betal, &
         ptcle_opt_betal_store, &
         ptcle_lamb_betal_store,&
         ptcle_nlamb_betal_store,&
         
         
         ptcle_material, &
         ptcle_param_mie3, &
         ptcle_param_hm,&
         ptcle_nphot_hm,&
         ptcle_param_baum,&
         ptcle_ilambda, &
         ptcle_ireflamb, &
         ptcle_Cext_reflamb,&
         ptcle_opt,&
         ptcle_nang_phasemat,&
         ptcle_u_phasemat,&
         ptcle_phasemat,&
         do_rt, &
         nang_max_store, &
         nlamb_max_store, &
         ptcle_opt_store, &
         ptcle_nlamb_phasemat_store, &
         ptcle_nang_phasemat_store, &
         ptcle_lamb_phasemat_store, &
         ptcle_u_phasemat_store, &
         ptcle_phasemat_store, &
         flag_store_opt, &
         nang_max_new, &
         nlamb_max_new, &
         ptcle_flag_new_opt, &
         ptcle_opt_new, &
         ptcle_nlamb_phasemat_new, &
         ptcle_nang_phasemat_new, &
         ptcle_lamb_phasemat_new, &
         ptcle_u_phasemat_new, &
         ptcle_phasemat_new, &
         flag_new_opt, &
         nang_max, &
         flag_baum,&
         flag_mie, &
         flag_hm, &
         nmu_baum

    USE MUTILITY, ONLY : LININTPOL
    !USE nr, ONLY : gauleg
    USE nm_nm, ONLY : nm_gaussNodes

    IMPLICIT NONE

    INTEGER :: i, j, k

    INTEGER :: icase
    INTEGER :: nlamb

    INTEGER, ALLOCATABLE             :: ilamb_ptcle(:)   ! to keep track of the lambda index
    !                                                      for a given ptcle

    INTEGER :: ncase
    INTEGER, ALLOCATABLE             :: itype(:)     ! to keep track of the type index 
    INTEGER, ALLOCATABLE             :: ilambda(:)   ! to keep track of the lambda index
    CHARACTER(LEN=max_len), ALLOCATABLE :: case_material(:)
    CHARACTER(LEN=max_len), ALLOCATABLE :: case_opt_model(:)
    CHARACTER(LEN=max_len), ALLOCATABLE :: case_ptcle_type(:)
    COMPLEX(KIND=DP), ALLOCATABLE :: case_refind(:)
    REAL(KIND=DP), ALLOCATABLE    :: case_lambda(:)
    REAL(KIND=DP), ALLOCATABLE    :: case_param_mie3(:,:)
    REAL(KIND=DP), ALLOCATABLE    :: case_param_hm(:,:)
    INTEGER, ALLOCATABLE          :: case_nphot_hm(:)
    REAL(KIND=DP), ALLOCATABLE    :: case_param_baum(:)
    INTEGER, ALLOCATABLE          :: case_nang_phasemat(:)
    REAL(KIND=DP), ALLOCATABLE    :: case_opt(:,:)
    REAL(KIND=DP), ALLOCATABLE    :: case_phasemat(:, :, :)
    REAL(KIND=DP), ALLOCATABLE    :: case_u_phasemat(:, :)
    INTEGER, ALLOCATABLE          :: case_nelem(:)

    !INTEGER, PARAMETER :: get_mie3_nangle = 901 
    !INTEGER, PARAMETER :: get_mie3_nangle = 181
    !INTEGER, PARAMETER :: get_hm_nangle = 901 

    !    INTEGER, PARAMETER :: get_mie3_nangle = 2800
    INTEGER, PARAMETER :: get_mie3_nangle = 2800
    INTEGER, PARAMETER :: get_hm_nangle = 901 

    INTEGER :: nang_tmp
    INTEGER :: nlamb_tmp
    REAL(kind=dp), ALLOCATABLE :: phasemat_tmp(:,:,:)
    REAL(kind=dp), ALLOCATABLE :: mu_tmp(:)
    REAL(kind=dp), ALLOCATABLE :: w_tmp(:)

    !------------

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute new optical properties
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get nlamb_max_new
    !     nang_max_new
    !     ncase : number of case for which to actually compute 
    !             new optical properties over all types and all wavelengths
    ALLOCATE(ptcle_flag_new_opt(nptcle))
    ptcle_flag_new_opt(:) = .FALSE.
    flag_new_opt =.FALSE.
    ncase = 0 

    nang_max_new = 0
    if (flag_mie) nang_max_new = get_mie3_nangle
    if ((flag_hm).and.(nang_max_new.lt.get_hm_nangle)) nang_max_new = get_hm_nangle
    if ((flag_baum).and.(nang_max_new.lt.nmu_baum)) nang_max_new = nmu_baum

    do i = 1, nptcle
       nlamb = 0
       !-- Line added on 19/02/14 
       if (ptcle_use_betal(i).eqv..false.) then
          if (ptcle_opt_interp(i).eqv..false.) then
             DO j = 1, nlambda
                IF (ptcle_opt_exist(i,j) .eqv. .FALSE.) then
                   nlamb = nlamb + 1
                   ncase = ncase + 1
                   ptcle_flag_new_opt(i) = .TRUE.
                endif
             ENDDO
             IF ((do_rt).and.(ireflamb.eq.0)) THEN
                IF (ptcle_opt_exist_reflamb(i) .eqv. .FALSE.) then
                   nlamb = nlamb + 1
                   ncase = ncase + 1
                   ptcle_flag_new_opt(i) = .TRUE.
                ENDIF
             ENDIF
             if (nlamb_max_new.lt.nlamb) nlamb_max_new = nlamb  
          end if  !!!!ptcle_opt_interp
       endif  ! if no betal read

    enddo
    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. get_opt) ',ncase,' new optical prop. will be computed'   
    ENDIF
    if (ncase.gt.0) flag_new_opt =.TRUE.

    if (flag_new_opt) THEN

       ! MCOMMON variables
       ALLOCATE(ptcle_opt_new(nptcle, ptcle_nopt, nlamb_max_new), &
            ptcle_nlamb_phasemat_new(nptcle), & 
            ptcle_nang_phasemat_new(nptcle,nlamb_max_new), &
            ptcle_lamb_phasemat_new(nptcle, nlamb_max_new), &
            ptcle_u_phasemat_new(nptcle, nang_max_new, nlamb_max_new), &
            ptcle_phasemat_new(nptcle, 6, nang_max_new, nlamb_max_new))
       ptcle_opt_new(:,:,:)          = undef_dp
       ptcle_nlamb_phasemat_new(:)   = undef_i
       ptcle_nang_phasemat_new(:,:)  = undef_i
       ptcle_lamb_phasemat_new(:, :) = undef_dp
       ptcle_u_phasemat_new(:,:,:)   = undef_dp
       ptcle_phasemat_new(:,:,:,:)   = undef_dp

       ! create arrays that contain all cases (no matter the type or wavelengths) 
       ! in order to have a single loop to parallelize
       ! that will contains ONLY cases for which we actually have to compute
       ! optical properties 
       ALLOCATE(itype(ncase), &
            ilambda(ncase), &
            case_material(ncase), &
            case_lambda(ncase), &
            case_opt_model(ncase), &
            case_ptcle_type(ncase),  &
            case_refind(ncase), &
            case_param_mie3(ncase,mie3_npar), &
            case_param_hm(ncase,hm_nmaxpar),&
            case_nphot_hm(ncase),&
            case_param_baum(ncase),&
            case_nang_phasemat(ncase), &
            case_opt(ncase,ptcle_nopt), &
            case_phasemat(ncase, 6, nang_max_new), &
            case_u_phasemat(ncase, nang_max_new),&
            case_nelem(ncase))
       ! init arrays
       itype(:)               = -1
       ilambda(:)             = -1
       case_material(:)       = undef_c
       case_lambda(:)         = undef_dp
       case_opt_model(:)      = undef_c
       case_ptcle_type(:)     = undef_c
       case_refind(:)         = CMPLX(undef_dp,undef_dp)
       case_param_mie3(:,:)   = undef_dp
       case_param_baum(:)     = undef_dp
       case_param_hm(:,:)     = undef_dp
       case_nphot_hm(:)       = undef_i
       case_nang_phasemat(:)  = undef_i
       case_opt(:,:)          = undef_dp
       case_phasemat(:,:,:)   = undef_dp
       case_u_phasemat(:,:)   = undef_dp
       case_nelem(:)          = undef_i 
       ! implement input "case" arrays
       icase = 0  
       DO i = 1, nptcle
          !---- line added on 19/02/14
          if (ptcle_use_betal(i).eqv..false.) then
             if (ptcle_opt_interp(i).eqv..false.) then
                IF  ((do_rt).and.(ireflamb.eq.0)) THEN
                   IF (ptcle_opt_exist_reflamb(i) .eqv. .FALSE.) THEN
                      icase = icase + 1
                      itype(icase)                   = i
                      ilambda(icase)                 = 0 ! index 0 for lambda means reflamb 
                      case_ptcle_type(icase)         = ptcle_type(i) 
                      case_material(icase)           = ptcle_material(i)
                      case_opt_model(icase)          = ptcle_opt_model(i)
                      case_lambda(icase)             = reflamb
                      case_param_mie3(icase,:)       = ptcle_param_mie3(i,:)
                      case_param_hm(icase,:)         = ptcle_param_hm(i,:)
                      case_nphot_hm(icase)           = ptcle_nphot_hm(i)
                      case_param_baum(icase)         = ptcle_param_baum(i)
                   ENDIF
                ENDIF
                DO j = 1, nlambda
                   IF (ptcle_opt_exist(i,j) .eqv. .FALSE.) THEN
                      icase = icase + 1
                      itype(icase)               = i
                      ilambda(icase)             = j
                      case_ptcle_type(icase)     = ptcle_type(i) 
                      case_material(icase)       = ptcle_material(i)
                      case_opt_model(icase)      = ptcle_opt_model(i)
                      case_lambda(icase)         = wlambda(j)
                      case_param_mie3(icase,:)   = ptcle_param_mie3(i,:)
                      case_param_hm(icase,:)     = ptcle_param_hm(i,:)
                      case_nphot_hm(icase)       = ptcle_nphot_hm(i)
                      case_param_baum(icase)     = ptcle_param_baum(i)
                   ENDIF
                ENDDO
             ENDIF  ! No ptcle_opt_interp
          ENDIF     ! No ptcle_use_betal

       ENDDO

       IF (verbose) WRITE (*,*) ''
       ! run optical properties computation
       DO icase = 1, ncase
          IF (verbose) THEN
             WRITE (*, FMT='(A49,A20,A2,F12.5,1x,A8)') & 
                  '  (sub. get_opt) compute optical properties for ', &
                  TRIM(ADJUSTL(case_ptcle_type(icase))), ' @',case_lambda(icase),' microns'             
          ENDIF
          SELECT CASE (case_opt_model(icase))
          case('mie3') 
             CALL GET_REFIND(case_material(icase), case_lambda(icase),  case_refind(icase))
             case_nang_phasemat(icase)  = get_mie3_nangle                  
             CALL CALL_MIE3(case_param_mie3(icase,:), case_nang_phasemat(icase),  case_refind(icase), case_lambda(icase), &
                  case_ptcle_type(icase), case_opt(icase,:), &
                  case_phasemat(icase, :, 1:case_nang_phasemat(icase)), case_u_phasemat(icase, 1:case_nang_phasemat(icase)) )
             case_nelem(icase) = 4
          case('baum') 
             ! note : for this one, case_nang_phasemat(icase) is an output
             CALL GET_BAUM(case_param_baum(icase), case_lambda(icase), nang_max_new, case_opt(icase,:),   &
                  case_nang_phasemat(icase), case_phasemat(icase, :, :), case_u_phasemat(icase, :))
             case_nelem(icase) = 6
          case('phm')                   
             CALL GET_REFIND(case_material(icase), case_lambda(icase),  case_refind(icase))
             case_nang_phasemat(icase)  = get_hm_nangle   
             CALL CALL_PHM(case_nphot_hm(icase), case_param_hm(icase,:), case_nang_phasemat(icase),  &
                  case_refind(icase), case_lambda(icase), &
                  case_ptcle_type(icase), case_opt(icase,:), &
                  case_phasemat(icase, :, 1:case_nang_phasemat(icase)), &
                  case_u_phasemat(icase, 1:case_nang_phasemat(icase)) )
             case_nelem(icase) = 6
          case('rhm')                   
             CALL GET_REFIND(case_material(icase), case_lambda(icase),  case_refind(icase))
             case_nang_phasemat(icase)  = get_hm_nangle   
             CALL CALL_RHM(case_nphot_hm(icase),case_param_hm(icase,:), case_nang_phasemat(icase),  &
                  case_refind(icase), case_lambda(icase), &
                  case_ptcle_type(icase), case_opt(icase,:), &
                  case_phasemat(icase, :, 1:case_nang_phasemat(icase)), &
                  case_u_phasemat(icase, 1:case_nang_phasemat(icase)) )
             case_nelem(icase) = 6
          case('ihm')                   
             CALL GET_REFIND(case_material(icase), case_lambda(icase),  case_refind(icase))
             case_nang_phasemat(icase)  = get_hm_nangle   
             CALL CALL_IHM(case_nphot_hm(icase),case_param_hm(icase,:), case_nang_phasemat(icase),  &
                  case_refind(icase), case_lambda(icase), &
                  case_ptcle_type(icase), case_opt(icase,:), &
                  case_phasemat(icase, :, 1:case_nang_phasemat(icase)), &
                  case_u_phasemat(icase, 1:case_nang_phasemat(icase)) )
             case_nelem(icase) = 6
          END SELECT
       ENDDO ! loop on cases

       ! implement result in common arrays
       ALLOCATE (ilamb_ptcle(nptcle))
       ilamb_ptcle(:) = 0
       ptcle_nlamb_phasemat_new(:)  = 0
       ptcle_nang_phasemat_new(:,:) = 0
       DO icase = 1, ncase

          ilamb_ptcle(itype(icase)) = ilamb_ptcle(itype(icase)) + 1  

          ! keep track of the index in OPT_NEW to find wlambda()
          if ( (ptcle_opt_exist_reflamb(itype(icase)).eqv..FALSE.) .and. &
               (ilambda(icase).eq.ireflamb) ) then
             ptcle_ireflamb(itype(icase)) = ilamb_ptcle(itype(icase))
          endif
          if (ilambda(icase).ne.0) then
             ptcle_ilambda(itype(icase),ilambda(icase)) = ilamb_ptcle(itype(icase))
          endif

          ptcle_nelem(itype(icase)) = case_nelem(icase)

          ptcle_opt_new(itype(icase),:,ilamb_ptcle(itype(icase)))         = case_opt(icase,:)
          ptcle_nlamb_phasemat_new(itype(icase))                          = ptcle_nlamb_phasemat_new(itype(icase)) + 1
          ptcle_lamb_phasemat_new(itype(icase),ilamb_ptcle(itype(icase))) = case_lambda(icase)
          ptcle_nang_phasemat_new(itype(icase),ilamb_ptcle(itype(icase))) = case_nang_phasemat(icase)
          write(*,*) 'icase, ptcle_nang_phasemat_new', icase,&
               ptcle_nang_phasemat_new(itype(icase),ilamb_ptcle(itype(icase)))

          ptcle_u_phasemat_new(itype(icase),1:case_nang_phasemat(icase),ilamb_ptcle(itype(icase))) = &
               case_u_phasemat(icase, 1:case_nang_phasemat(icase))
          ptcle_phasemat_new(itype(icase),:,1:case_nang_phasemat(icase),ilamb_ptcle(itype(icase))) = &
               case_phasemat(icase, :, 1:case_nang_phasemat(icase))
       ENDDO ! loop on cases

       DEALLOCATE(itype, case_material, case_lambda, &
            case_opt_model, case_ptcle_type,  &
            case_refind, case_param_mie3, &
            case_param_hm, case_nphot_hm, &
            case_param_baum, &
            case_nang_phasemat, case_opt, case_phasemat, &
            case_u_phasemat, ilamb_ptcle)

    endif ! end compute new opt

    ! get nang_max and alllocate OPT arrays
    nang_max = 0
    if (flag_new_opt) then
       if (nang_max.lt.nang_max_new) nang_max = nang_max_new
    endif
    if (flag_store_opt) then
       if (nang_max.lt.nang_max_store) nang_max = nang_max_store
    endif

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  GET PTCLE OPTICAL PROPERTIES
    !  TO BE USED IN THE REST OF THE RUN (GET_BETAL, DO_RT ...)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (OPT_ONLY.EQV..FALSE.) THEN

       IF (verbose) THEN
          WRITE (*,*) ''
          WRITE (*, *)'  (sub. get_opt) Set-up the optical properties variables for Betal expansion and RT'
       ENDIF


       ALLOCATE(ptcle_Cext_reflamb(nptcle),&
            ptcle_opt(nptcle, ptcle_nopt, nlambda), &
            ptcle_nang_phasemat(nptcle,nlambda) , &
            ptcle_u_phasemat(nptcle, nang_max, nlambda),&
            ptcle_phasemat(nptcle, 6, nang_max, nlambda))

       ptcle_Cext_reflamb(:)    = undef_dp
       ptcle_opt(:,:,:)         = undef_dp
       ptcle_nang_phasemat(:,:) = undef_i
       ptcle_u_phasemat(:,:,:)  = undef_dp
       ptcle_phasemat(:,:,:,:)  = 0.0_dp


       ALLOCATE(mu_tmp(nang_max_store), w_tmp(nang_max_store), phasemat_tmp(6,nang_max_store,nlamb_max_store))

       DO i = 1, nptcle
          !---- line added on 19/02/14
          !----- line modified  on 21/02/14
          if (ptcle_use_betal(i).eqv..true.) then               
             ! get Cext for reference wavelengths
             if (ptcle_opt_exist_reflamb(i)) then
                ! if the wavelength is stored there is no need for interpolation
                ptcle_Cext_reflamb(i) = ptcle_opt_betal_store(i,1,ptcle_ireflamb(i))
             else               
                write(*,*) ''
                WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                WRITE(*,*)               ' (in sub. get_opt)  : ERROR '
                write(*,*)               '   For particle ', TRIM(ADJUSTL(ptcle_type(i)))
                write(*,FMT='(1x,A,F12.5,2x,A)') '   the working wavelenght '
                write(*,*)               '   must set the extinction coefficient at 550 as '
                write(*,*)               '   reference'
                WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*) ''
                STOP
             endif


             DO j = 1, nlambda

                ! here we test wheteher the required wavelength is inside the stored range
                if ( (wlambda(j).lt.MINVAL(ptcle_lamb_betal_store(i,1:ptcle_nlamb_betal_store(i)))) .or. &
                     (wlambda(j).gt.MAXVAL(ptcle_lamb_betal_store(i,1:ptcle_nlamb_betal_store(i))))  ) then                  

                   write(*,*) ''
                   WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   WRITE(*,*)               ' (in sub. get_opt)  : ERROR '
                   write(*,*)               '   For particle ', TRIM(ADJUSTL(ptcle_type(i)))
                   write(*,FMT='(1x,A,F12.5,2x,A)') '   the working wavelenght ', wlambda(j), 'microns'
                   write(*,*)               '   is outside the wavelength range of stored'
                   write(*,*)               '   optical properties'
                   WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   write(*,*) ''
                   STOP

                else                        

                   if (ptcle_opt_exist(i,j)) then
                      ! if the wavelength is stored there is no need for interpolation
                      ptcle_opt(i,:,j)  = ptcle_opt_betal_store(i,:,ptcle_ilambda(i,j))
                   else
                      write(*,*) ''
                      WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      WRITE(*,*)               ' (in sub. get_opt)  : ERROR '
                      write(*,*)               '   For particle ', TRIM(ADJUSTL(ptcle_type(i)))
                      write(*,FMT='(1x,A,F12.5,2x,A)') '   the working wavelenght ', wlambda(j), 'microns'
                      write(*,*)               '   require optical properties such as Cext, albedo '
                      write(*,*)               '   simple difusion, ...'
                      WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      write(*,*) ''
                      STOP

                   endif  !ptcle_opt_exist(i,j)

                endif  !the required wavelength is inside the stored range
             END DO !nlambda



          else
             if (ptcle_opt_interp(i)) then

                !------------------------------------
                ! INTERPOLATION OF OPTICAL PROPERTIES

                if (DO_RT) then
                   ! get Cext for reference wavelengths
                   if (ptcle_opt_exist_reflamb(i)) then
                      ! if the wavelength is stored there is no need for interpolation
                      ptcle_Cext_reflamb(i) = ptcle_opt_store(i,1,ptcle_ireflamb(i))
                   else                   


                      ptcle_Cext_reflamb(i) =  LININTPOL(ptcle_opt_store(i,1,1:ptcle_nlamb_phasemat_store(i)),&
                           ptcle_lamb_phasemat_store(i,1:ptcle_nlamb_phasemat_store(i)), &
                           ptcle_nlamb_phasemat_store(i), &
                           reflamb)


                   endif
                endif

                ! get Cext, SSA and g 
                DO j = 1, nlambda

                   ! here we test wheteher the required wavelength is inside the stored range
                   if ( (wlambda(j).lt.MINVAL(ptcle_lamb_phasemat_store(i,1:ptcle_nlamb_phasemat_store(i)))) .or. &
                        (wlambda(j).gt.MAXVAL(ptcle_lamb_phasemat_store(i,1:ptcle_nlamb_phasemat_store(i))))  ) then                   

                      write(*,*) ''
                      WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      WRITE(*,*)               ' (in sub. get_opt)  : ERROR '
                      write(*,*)               '   For particle ', TRIM(ADJUSTL(ptcle_type(i)))
                      write(*,FMT='(1x,A,F12.5,2x,A)') '   the working wavelenght ', wlambda(j), 'microns'
                      write(*,*)               '   is outside the wavelength range of stored'
                      write(*,*)               '   optical properties'
                      WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      write(*,*) ''
                      STOP

                   else

                      if (ptcle_opt_exist(i,j)) then
                         ! if the wavelength is stored there is no need for interpolation
                         ptcle_opt(i,:,j)  = ptcle_opt_store(i,:,ptcle_ilambda(i,j))
                      else
                         ! Cext, SSA, g for all wavelengths
                         DO k = 1, 3
                            ptcle_opt(i,k,j) =  LININTPOL(ptcle_opt_store(i,k,1:ptcle_nlamb_phasemat_store(i)),&
                                 ptcle_lamb_phasemat_store(i,1:ptcle_nlamb_phasemat_store(i)), &
                                 ptcle_nlamb_phasemat_store(i), &
                                 wlambda(j))
                         END DO
                      endif

                   endif

                END DO

                ! Set the phase matrix on a common mu grid over the different wavelengths
                ! and interpolate 
                nlamb_tmp = ptcle_nlamb_phasemat_store(i)
                nang_tmp  = MAXVAL(ptcle_nang_phasemat_store(i, 1:ptcle_nlamb_phasemat_store(i)))
                !CALL gauleg(-1.0_dp, +1.0_dp, mu_tmp(2:nang_tmp-1), w_tmp(2:nang_tmp-1))
                !write(*,*) , 'artdeco_get_opt (  gauleg )   before mu_tmp , w_tmp ', mu_tmp(2:nang_tmp-1),  w_tmp(2:nang_tmp-1)

                CALL nm_gaussNodes(-1.0_dp, +1.0_dp, mu_tmp(2:nang_tmp-1), w_tmp(2:nang_tmp-1))
                !write(*,*) , 'artdeco_get_opt (  gauleg )   after mu_tmp ,nang_tmp  ', mu_tmp(2:nang_tmp-1),nang_tmp

                ! set the extremum values by hand cause gauss grid does not include
                mu_tmp(1)        = -1.0_dp
                mu_tmp(nang_tmp) = +1.0_dp

                if (ptcle_hg(i)) then

                   ! H-G phase function
                   do j = 1, nlambda 
                      ptcle_nang_phasemat(i,j) = nang_tmp
                      ptcle_u_phasemat(i, 1:nang_tmp, j) = mu_tmp(1:nang_tmp)
                      do k = 1, ptcle_nang_phasemat(i,j) 
                         ptcle_phasemat(i, 1, k, j) = HG(ptcle_u_phasemat(i,k,j), ptcle_opt(i,3,j))
                      end do
                   end do

                else

                   ! We first need to interpolate on the common mu grid
                   if (nmat.ge.1) then
                      ! F11
                      DO j = 1, nlamb_tmp
                         DO k = 1, nang_tmp
                            phasemat_tmp(1,k,j) = LININTPOL(ptcle_phasemat_store(i, 1, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                         END DO
                      END DO
                   end if
                   if (nmat.ge.3) then
                      ! F22, F33, F21
                      DO j = 1, nlamb_tmp
                         DO k = 1, nang_tmp
                            phasemat_tmp(2,k,j) = LININTPOL(ptcle_phasemat_store(i, 2, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                            phasemat_tmp(3,k,j) = LININTPOL(ptcle_phasemat_store(i, 3, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                            phasemat_tmp(5,k,j) = LININTPOL(ptcle_phasemat_store(i, 5, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                         END DO
                      END DO
                   end if
                   if (nmat.eq.4) then
                      ! F44 & F34
                      DO j = 1, nlamb_tmp
                         DO k = 1, nang_tmp
                            phasemat_tmp(4,k,j) = LININTPOL(ptcle_phasemat_store(i, 4, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                            phasemat_tmp(6,k,j) = LININTPOL(ptcle_phasemat_store(i, 6, 1:ptcle_nang_phasemat_store(i,j), j),&
                                 ptcle_u_phasemat_store(i, 1:ptcle_nang_phasemat_store(i,j), j), &
                                 ptcle_nang_phasemat_store(i,j), &
                                 mu_tmp(k))
                         END DO
                      END DO
                   end if

                   ! now we can interpolate the phase matrix on working lambda
                   DO j = 1, nlambda

                      if (ptcle_opt_exist(i,j)) then

                         ! if the wavelength is stored there is no need for interpolation
                         ptcle_nang_phasemat(i,j) = ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j))
                         ptcle_u_phasemat(i,1:ptcle_nang_phasemat(i,j),j)  = &
                              ptcle_u_phasemat_store(i,1:ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))
                         ptcle_phasemat(i, :, 1:ptcle_nang_phasemat(i,j), j) = &
                              ptcle_phasemat_store(i, :, 1:ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))

                      else

                         ptcle_nang_phasemat(i,j) = nang_tmp

                         ptcle_u_phasemat(i, 1:nang_tmp, j) = mu_tmp(1:nang_tmp)

                         if (nmat.ge.1) then
                            ! F11
                            DO k = 1, nang_tmp
                               ptcle_phasemat(i, 1, k, j) = LININTPOL( phasemat_tmp(1,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                            END DO
                         end if
                         if (nmat.ge.3) then
                            ! F11
                            DO k = 1, nang_tmp
                               ptcle_phasemat(i, 2, k, j) = LININTPOL( phasemat_tmp(2,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                               ptcle_phasemat(i, 3, k, j) = LININTPOL( phasemat_tmp(3,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                               ptcle_phasemat(i, 5, k, j) = LININTPOL( phasemat_tmp(5,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                            END DO
                         end if
                         if (nmat.eq.4) then
                            ! F11
                            DO k = 1, nang_tmp
                               ptcle_phasemat(i, 4, k, j) = LININTPOL( phasemat_tmp(4,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                               ptcle_phasemat(i, 6, k, j) = LININTPOL( phasemat_tmp(6,k,1:nlamb_tmp), &
                                    ptcle_lamb_phasemat_store(i, 1:nlamb_tmp), &
                                    nlamb_tmp, &
                                    wlambda(j))
                            END DO
                         end if

                      endif

                   ENDDO! lambda loop

                endif! H-G

             else

                !---------------------------------------
                ! NO INTERPOLATION OF OPTICAL PROPERTIES

                if (DO_RT) then
                   if (ptcle_opt_exist_reflamb(i)) then
                      ptcle_Cext_reflamb(i) = ptcle_opt_store(i,1,ptcle_ireflamb(i))
                   else
                      ptcle_Cext_reflamb(i) = ptcle_opt_new(i,1,ptcle_ireflamb(i))
                   endif
                endif

                DO j = 1, nlambda
                   if (ptcle_opt_exist(i,j)) then
                      ptcle_opt(i,:,j)         = ptcle_opt_store(i,:,ptcle_ilambda(i,j))
                      ptcle_nang_phasemat(i,j) = ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j))

                      ptcle_u_phasemat(i,1:ptcle_nang_phasemat(i,j),j)  = &
                           ptcle_u_phasemat_store(i,1:ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))
                   else
                      ptcle_opt(i,:,j)         = ptcle_opt_new(i,:,ptcle_ilambda(i,j))
                      ptcle_nang_phasemat(i,j) = ptcle_nang_phasemat_new(i,ptcle_ilambda(i,j))
                      !write(*,*) 'ptcle_ilambda,ptcle_nang_phasemat_store,  ptcle_u_phasemat_store',ptcle_ilambda, &
                      !ptcle_nang_phasemat_new(i,ptcle_ilambda(i,j)), &
                      !ptcle_u_phasemat_new(i,1:ptcle_nang_phasemat_new(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))

                      ptcle_u_phasemat(i,1:ptcle_nang_phasemat(i,j),j) = &
                           ptcle_u_phasemat_new(i,1:ptcle_nang_phasemat_new(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))
                      !write(*,*) 'i,ptcle_nang_phasemat, ptcle_u_phasemat', i, ptcle_nang_phasemat(i,j),&
                      !ptcle_u_phasemat(i,1:ptcle_nang_phasemat(i,j),j)
                   endif
                ENDDO

                ! The phase matrix
                if (ptcle_hg(i)) then
                   ! H-G phase function
                   do j = 1, nlambda 
                      do k = 1, ptcle_nang_phasemat(i,j) 
                         ptcle_phasemat(i, 1, k, j) = HG(ptcle_u_phasemat(i,k,j), ptcle_opt(i,3,j))
                      end do
                   end do
                else
                   DO j = 1, nlambda
                      if (ptcle_opt_exist(i,j)) then
                         ptcle_phasemat(i, :, 1:ptcle_nang_phasemat(i,j), j) = &
                              ptcle_phasemat_store(i, :, 1:ptcle_nang_phasemat_store(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))
                      else
                         ptcle_phasemat(i, :, 1:ptcle_nang_phasemat(i,j), j) = &
                              ptcle_phasemat_new(i, :, 1:ptcle_nang_phasemat_new(i,ptcle_ilambda(i,j)), ptcle_ilambda(i,j))
                      endif
                   END DO
                endif

             endif ! interp or not
          endif !No read betal
       end DO ! loop on ptcle

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! DO i = 1, nptcle
       !    DO j = 1, nlambda 
       !       print*, ' get_opt'
       !       print*, TRIM(ADJUSTL(ptcle_type(i))),'  ',lambda(j)
       !       DO k = 1, ptcle_nang_phasemat(i,j)
       !          print*, ptcle_u_phasemat(i,k,j), &
       !               ptcle_phasemat(i, 1, k, j), &
       !               ptcle_phasemat(i, 4, k, j),  &
       !               ptcle_phasemat(i, 5, k, j),  &
       !               ptcle_phasemat(i, 6, k, j)
       !       ENDDO
       !    ENDDO
       ! END DO
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       DEALLOCATE(mu_tmp, w_tmp, phasemat_tmp)

    ENDIF

  END SUBROUTINE GET_OPT

  !----------------------------------------------------------------

  FUNCTION HG(mu, g)

    USE MCONSTANTS, only : dp

    IMPLICIT NONE

    REAL(kind=dp) :: HG
    REAL(kind=dp), INTENT(IN) :: mu
    REAL(kind=dp), INTENT(IN) :: g

    HG = ( 1.0_dp - g**2.0_dp ) / &
         (( 1.0_dp + g**2.0_dp - 2.0_dp*g*mu)**(3.0_dp/2.0_dp)) 

  END FUNCTION HG

  !----------------------------------------------------------------

  SUBROUTINE GET_REFIND(material, lamb, refind)

    USE MCOMMON, only : nmaterial, &
         material_name, &
         material_nlamb_refind, &
         material_refind
    USE MCONSTANTS, only : dp, max_len, undef_dp
    USE MUTILITY, only : LININTPOL

    IMPLICIT NONE

    CHARACTER(LEN=max_len), INTENT(IN) :: material
    REAL(kind=dp), INTENT(IN)       :: lamb
    COMPLEX(KIND=dp), INTENT(OUT)   :: refind

    ! local variable to call refwat and refice
    INTEGER       :: imaterial
    REAL(kind=dp) :: refind_re
    REAL(kind=dp) :: refind_im
    INTEGER :: nlamb

    !-------
    refind_re = undef_dp
    refind_im = undef_dp

    imaterial = 1
    do while (TRIM(ADJUSTL(material_name(imaterial))) .ne. TRIM(ADJUSTL(material)))
       imaterial = imaterial + 1
    end do

    !write(*,*) TRIM(ADJUSTL(material_name(imaterial))), TRIM(ADJUSTL(material))

    nlamb = material_nlamb_refind(imaterial)
    if ( (lamb .lt. minval(material_refind(imaterial,1,1:nlamb))) .or.   &
         (lamb .gt. maxval(material_refind(imaterial,1,1:nlamb))) ) then
       WRITE(*,*) ''
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_refind)  : ERROR '
       WRITE(*,*) ' For material ', TRIM(ADJUSTL(material))
       WRITE(*,*) ' wavelengths ', lamb
       WRITE(*,*) ' is out of defined ref. ind. range' 
       WRITE(*,*) ' lambmin = ', minval(material_refind(imaterial,1,1:nlamb))
       WRITE(*,*) ' lambmax = ', maxval(material_refind(imaterial,1,1:nlamb))
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP                       
    endif

    ! note : interpolation are done on LOG(wavelengths)
    refind_re = LININTPOL( &
         material_refind(imaterial,2,1:nlamb), &        ! real. part.
         LOG10(material_refind(imaterial,1,1:nlamb)), & ! wavelengths
         nlamb, &                                       ! nwavelengths
         LOG10(lamb))                                   ! wavelength to get real. part for

    refind_im = LININTPOL( &
         material_refind(imaterial,3,1:nlamb), &        ! im. part.
         LOG10(material_refind(imaterial,1,1:nlamb)), & ! wavelengths
         nlamb, &                                       ! nwavelengths
         LOG10(lamb))                                   ! wavelength to get im. part for

    refind = CMPLX(refind_re, -refind_im)
    !write(*,*)refind

  END SUBROUTINE GET_REFIND

  !----------------------------------------------------------------

  SUBROUTINE GET_BAUM(param_baum, lamb, nmaxang, opt, nang, phasemat, u_phasemat)

    USE MCONSTANTS, only : dp, &
         ptcle_nopt
    USE MCOMMON, only : nlamb_baum, &
         nde_baum,  &
         nmu_baum,  &
         lamb_baum, & ! lamb_baum(nlamb_baum)
         de_baum,   & ! de_baum(nde_baum)
         mu_baum,   & ! mu_baum(nmu_baum)
         Cext_baum, & ! Cext_baum(nlamb_baum, nde_baum)
         ssa_baum,  & ! ssa_baum(nlamb_baum, nde_baum)
         g_baum,    & ! g_baum(nlamb_baum, nde_baum)
         ptcle_phasemat_baum  ! ptcle_phasemat_baum(nlamb_baum, nde_baum, 6, nmu_baum) : the phase matrix
    USE MUTILITY, only : BILININTPOL

    IMPLICIT NONE

    REAL(kind=dp), INTENT(IN)    :: param_baum
    REAL(kind=dp), INTENT(IN)    :: lamb
    INTEGER, INTENT(IN)          :: nmaxang
    REAL(kind=dp), INTENT(OUT)   :: opt(ptcle_nopt)
    INTEGER, INTENT(OUT)         :: nang
    REAL(kind=dp), INTENT(OUT)   :: phasemat(6,nmaxang)
    REAL(kind=dp), INTENT(OUT)   :: u_phasemat(nmaxang)

    integer :: i

    !-------

    opt(1) = BILININTPOL(Cext_baum, lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum) ! Cext
    opt(2) = BILININTPOL(ssa_baum, lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)  ! SSA
    opt(3) = BILININTPOL(g_baum, lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)    ! g : asymetry fcator
    nang   = nmu_baum

    if (nmaxang.lt.nang) then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. GET_BAUM) :  ERROR '
       WRITE(*,*) '  nmaxang should be greater or equal to nang' 
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
       STOP
    endif

    do i = 1, nang
       u_phasemat(i) = mu_baum(i)
       phasemat(1,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 1, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
       phasemat(2,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 2, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
       phasemat(3,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 3, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
       phasemat(4,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 4, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
       phasemat(5,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 5, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
       phasemat(6,i) = BILININTPOL(ptcle_phasemat_baum(:, :, 6, i), lamb_baum, de_baum, nlamb_baum, nde_baum, lamb, param_baum)
    end do

  END SUBROUTINE GET_BAUM

  !----------------------------------------------------------------

  SUBROUTINE CALL_MIE3(param_mie3, nang, refind, lamb, ptype, opt, phasemat, u_phasemat)

    USE MCONSTANTS, only : dp, &
         max_len, &
         mie3_npar, &
         ptcle_nopt

    USE MCOMMON, only :warning
    USE MUTILITY, only: PLINT

    IMPLICIT NONE

    REAL(kind=dp), INTENT(IN)    :: param_mie3(mie3_npar)
    INTEGER, INTENT(IN)          :: nang
    COMPLEX(kind=dp), INTENT(IN) :: refind
    REAL(kind=dp), INTENT(IN)    :: lamb
    CHARACTER(LEN=max_len), INTENT(IN) :: ptype
    REAL(kind=dp), INTENT(OUT)   :: opt(ptcle_nopt)
    REAL(kind=dp), INTENT(OUT)   :: phasemat(6,nang)
    REAL(kind=dp), INTENT(OUT)   :: u_phasemat(nang)

    ! local variable to call mie3
    INTEGER, PARAMETER :: mie3_nmiec = 10 ! number of element in the miec(:) array out of mie3
    ! input
    REAL(kind=dp) :: mie3_delta  
    REAL(kind=dp) :: mie3_cutoff
    INTEGER :: mie3_idis
    INTEGER :: mie3_nsubr
    INTEGER :: mie3_ngaur
    REAL(kind=dp) :: mie3_par1
    REAL(kind=dp) :: mie3_par2
    REAL(kind=dp) :: mie3_par3
    INTEGER :: mie3_igaussu
    INTEGER:: mie3_nangle
    !output
    REAL(kind=dp) :: mie3_F(4,nang)
    !                mie3_F(1,:) = F11 = F22
    !                mie3_F(2,:) = F21 = F12
    !                mie3_F(3,:) = F33 = F44
    !                mie3_F(4,:) = F34 = -F43
    REAL(kind=dp) :: mie3_miec(mie3_nmiec)
    ! mie3_miec(nptcle,mie3_nmiec)   : optical parameters output from mie3
    ! mie3_miec(1)  : Csca   : the average scattering cross section (microns^2)
    ! mie3_miec(2)  : Cext   : the average extinction cross section (microns^2)
    ! mie3_miec(3)  : Qsca   : the scattering efficiency factor 
    ! mie3_miec(4)  : Qext   : the extinction efficiency factor 
    ! mie3_miec(5)  : w0     : the single scattering albedo
    ! mie3_miec(6)  : G      : the average geometrical cross section (microns^2)
    ! mie3_miec(7)  : reff   : the average effective radius (microns)
    ! mie3_miec(8)  : xeff   : the average effective size parameter
    ! mie3_miec(9)  : numpar : the integrated number of particles
    ! mie3_miec(10) : volume : the average volume  (microns^3)
    REAL(kind=dp) :: mie3_u(nang)
    REAL(kind=dp) :: mie3_w(nang)

    REAL(kind=dp) :: res_plint

    !----------

    mie3_delta   = param_mie3(1)
    mie3_cutoff  = param_mie3(2)
    mie3_idis    = NINT(param_mie3(3))
    mie3_nsubr   = NINT(param_mie3(4))
    mie3_ngaur   = NINT(param_mie3(5))
    mie3_par1    = param_mie3(6)
    mie3_par2    = param_mie3(7)
    mie3_par3    = param_mie3(8)

    !mie3_igaussu = 0    ! if 1 : phase matrix will be computed over mie3_nangle Gauss points between 0 and 180. deg.
    !                     if 0 : phase matrix on nangle equidistant mie3_nangle between 0 and 180. deg.   
    mie3_igaussu = 1
    mie3_nangle  = nang

    IF (mie3_delta .lt. 0) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   In calling mie3 -delta- must be positive '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    ENDIF
    IF (mie3_cutoff .lt. 0) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   In calling mie3 -cutoff- must be positive '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    ENDIF
    if (REAL(refind) .le. 0.0_dp) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   at lamb = ', lamb
       WRITE(*,*) '   Real part of refractive index must be > 0'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    endif
    if (AIMAG(refind) .gt. 0.0_dp) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   at lamb = ', lamb
       WRITE(*,*) '   Imaginary part of refractive index must be <= 0'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    endif
    if ((mie3_idis .lt. 0) .or. (mie3_idis .gt. 8)) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   In calling mie3 size distribution index must be 0 <= idis <= 8'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    endif
    IF (mie3_nsubr .le. 0) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   In calling mie3 the number of size interval (nsubr) must be greater than 0 '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    ENDIF
    IF (mie3_ngaur .le. 0) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. get_betal)  : ERROR '
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   In calling mie3 the number of Gauss points must be greater than 0 '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    ENDIF

    call mie3(mie3_nangle, mie3_delta, mie3_cutoff       &
         , lamb, refind, mie3_nsubr , mie3_ngaur, mie3_idis    &
         , mie3_par1, mie3_par2, mie3_par3, mie3_igaussu, mie3_nangle, mie3_F, mie3_u, mie3_w, mie3_miec)

    if (mie3_miec(5) .gt. 1.0_dp) THEN
       if (warning)then
          WRITE(*,*) ' '
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (in sub. call_mie3)     WARNING'
          WRITE(*,*) '   For particle : ', ptype
          WRITE(*,*) '   at lamb = ', lamb
          WRITE(*,*) '   Scattering albedo can not be > 1.0 '
          WRITE(*,*) '      w0 = ', mie3_miec(5)
          WRITE(*,*) '   We reset it to 1.0'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' '
       endif
       mie3_miec(5) = 1.0_dp
    endif

    opt(1) =  mie3_miec(2)  ! Cext
    opt(2) =  mie3_miec(5)  ! SSA

    ! compute asymetry factor
    CALL PLINT( mie3_u(1:nang)*mie3_F(1,1:nang) / 2.0_dp, mie3_u(1:nang), nang, -1.0_dp, 1.0_dp, res_plint)
    opt(3) = res_plint  ! g : asymetry factor
    !opt(3) =  -32768!res_plint  ! g : asymetry factor

    phasemat(1, :) = mie3_F(1,:) ! F11
    phasemat(2, :) = mie3_F(1,:) ! F22
    phasemat(3, :) = mie3_F(3,:) ! F33
    phasemat(4, :) = mie3_F(3,:) ! F44
    phasemat(5, :) = mie3_F(2,:) ! F21
    phasemat(6, :) = mie3_F(4,:) ! F34
    u_phasemat(:)  = mie3_u

  END SUBROUTINE CALL_MIE3

  !----------------------------------------------------------------

  SUBROUTINE CALL_PHM(nphot, param_hm, nang, refind, lamb, ptype, opt, phasemat, u_phasemat)

    USE MCOMMON, only : warning
    USE MCONSTANTS, only : dp, &
         max_len, &
         ptcle_nopt, &
         deg2rad, &
         hm_nmaxpar, &
         min_size_param_hm, &
         xpi
    !USE nr, ONLY : SORT2
    USE nm_nm, ONLY : NM_SORT2


    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nphot
    REAL(kind=dp), INTENT(IN)    :: param_hm(hm_nmaxpar)
    INTEGER, INTENT(IN)          :: nang
    COMPLEX(kind=dp), INTENT(IN) :: refind
    REAL(kind=dp), INTENT(IN)    :: lamb
    CHARACTER(LEN=max_len), INTENT(IN) :: ptype
    REAL(kind=dp), INTENT(OUT)   :: opt(ptcle_nopt)
    REAL(kind=dp), INTENT(OUT)   :: phasemat(6,nang)
    REAL(kind=dp), INTENT(OUT)   :: u_phasemat(nang)

    ! local variables
    INTEGER, PARAMETER :: hm_nout = 8   ! number of values in hm_PHM_RES(:) 
    ! HM input
    REAL(kind=dp) :: hm_rm
    REAL(kind=dp) :: hm_im
    REAL(kind=dp) :: hm_ll
    REAL(kind=dp) :: hm_rr
    INTEGER       :: hm_nangle
    ! HM output
    REAL(kind=dp) :: hm_PHM_RES(hm_nout)
    !                hm_phm_res(1) = Number of scattered photons
    !                hm_phm_res(2) = Incident energy flux
    !                hm_phm_res(3) = scattered energy flux
    !                hm_phm_res(4) = diffracted energy flux
    !                hm_phm_res(5) = extinction cross section
    !                hm_phm_res(6) = scattering cross section
    !                hm_phm_res(7) = ssa
    !                hm_phm_res(8) = asymetry factor
    REAL(kind=dp), allocatable :: hm_PGT(:,:,:)
    REAL(kind=dp), allocatable :: hm_thetaout(:)

    REAL(kind=dp), allocatable :: tmp_u(:)
    REAL(kind=dp) :: size_param
    REAL(kind=dp) :: r_eq
    REAL(kind=dp) :: v_crystal

    INTEGER :: i

    !========

    hm_ll = param_hm(1)
    hm_rr = param_hm(2)
    ! test ray tracing method validity
    ! compute the radius of the equivalent volume sphere
    v_crystal = (3.0_dp*SQRT(3.0_dp)/2.0_dp * hm_rr**2.0_dp * hm_ll)
    r_eq      = (3.0_dp/4.0_dp*v_crystal/xpi)**(1.0_dp/3.0_dp)
    size_param = 2.0_dp * xpi * r_eq / lamb
    if (size_param .lt. min_size_param_hm) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. CALL_PHM)  : ERROR '
       WRITE(*,*) '   For the ray tracing method to be valid, '
       WRITE(*,*) '   the size parameter must be greater than', min_size_param_hm
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   at lamb = ', lamb
       WRITE(*,*) '   2*PI*r_eq/lamb = ', size_param
       WRITE(*,*) '   r_eq= ', r_eq
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP
    endif

    hm_rm = REAL(refind)
    hm_im = AIMAG(refind)
    hm_nangle = nang
    ALLOCATE(hm_PGT(4,4,hm_nangle),  &
         hm_thetaout(hm_nangle), &
         tmp_u(hm_nangle) )

    CALL hm_phm2p4(nphot, lamb, hm_rm, hm_im,  &
         hm_ll, hm_rr, hm_nangle,         & 
         hm_PHM_RES, hm_thetaout, hm_PGT)   

    if (hm_PHM_RES(7) .gt. 1.0_dp) THEN
       if (warning) then
          WRITE(*,*) ' '
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (in sub. call_PHM)     WARNING'
          WRITE(*,*) '   Scattering albedo can not be > 1.0 '
          WRITE(*,*) '   For particles ', ptype
          WRITE(*,*) '      w0 = ',hm_PHM_RES(7)
          WRITE(*,*) '   We reset it to 1.0'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' '
       endif
       hm_PHM_RES(7) = 1.0_dp
    endif

!!!! line modified in 27/03/2015 for sort using another algo
    ! sort in increasing u
    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(1,1,:))  ! F11
    call NM_SORT2(tmp_u, hm_PGT(1,1,:))  ! F11

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,2,:))  ! F22
    call NM_SORT2(tmp_u, hm_PGT(2,2,:))  ! F22

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,3,:))  ! F33
    call NM_SORT2(tmp_u, hm_PGT(3,3,:))  ! F33

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(4,4,:))  ! F44
    call NM_SORT2(tmp_u, hm_PGT(4,4,:))  ! F44

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,1,:))  ! F21
    call NM_SORT2(tmp_u, hm_PGT(2,1,:))  ! F21

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,4,:))  ! F34
    call NM_SORT2(tmp_u, hm_PGT(3,4,:))  ! F34


    ! write phase matrix in the common variable
    DO i= 1, hm_nangle
       u_phasemat(i)  = tmp_u(i)
       phasemat(1, i) = hm_PGT(1,1,i)               ! F11
       phasemat(2, i) = hm_PGT(2,2,i)*hm_PGT(1,1,i) ! F22
       phasemat(3, i) = hm_PGT(3,3,i)*hm_PGT(1,1,i) ! F33
       phasemat(4, i) = hm_PGT(4,4,i)*hm_PGT(1,1,i) ! F44
       phasemat(5, i) = hm_PGT(2,1,i)*hm_PGT(1,1,i) ! F21
       phasemat(6, i) = hm_PGT(3,4,i)*hm_PGT(1,1,i) ! F34
    ENDDO

    ! The computation is done over angle domains [thetamin:thetamax]
    ! so that we do not have 0 and 180 deg points in the final grid.
    ! We force the last and first angle to be
    ! 0 and 180 deg. 
    if(warning)then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_PHM)     WARNING'
       WRITE(*,*) ' The computation is done over angle domains [thetamin:thetamax]'
       WRITE(*,*) ' so that we do not have 0 and 180 deg points in the final grid.'
       WRITE(*,*) ' We force the last and first angle to be'
       WRITE(*,*) ' 0 and 180 deg. '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
    endif
    u_phasemat(1)         = -1.0_dp
    u_phasemat(hm_nangle) = 1.0_dp

    opt(1) = hm_phm_res(5) ! Cext
    opt(2) = hm_phm_res(7) ! ssa
    opt(3) = hm_phm_res(8) ! asymetry factor

    DEALLOCATE(hm_PGT, hm_thetaout, tmp_u)

  END SUBROUTINE CALL_PHM

  !----------------------------------------------------------------

  SUBROUTINE CALL_RHM(nphot, param_hm, nang, refind, lamb, ptype, opt, phasemat, u_phasemat)

    USE MCONSTANTS, only : dp, &
         max_len, &
         ptcle_nopt, &
         deg2rad, &
         hm_nmaxpar, &
         min_size_param_hm, &
         xpi
    USE MCOMMON, ONLY : warning
    !USE nr, ONLY : SORT2
    USE nm_nm, ONLY : NM_SORT2

    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nphot
    REAL(kind=dp), INTENT(IN)    :: param_hm(hm_nmaxpar)
    INTEGER, INTENT(IN)          :: nang
    COMPLEX(kind=dp), INTENT(IN) :: refind
    REAL(kind=dp), INTENT(IN)    :: lamb
    CHARACTER(LEN=max_len), INTENT(IN) :: ptype
    REAL(kind=dp), INTENT(OUT)   :: opt(ptcle_nopt)
    REAL(kind=dp), INTENT(OUT)   :: phasemat(6,nang)
    REAL(kind=dp), INTENT(OUT)   :: u_phasemat(nang)

    ! local variables
    INTEGER, PARAMETER :: hm_nout = 8   ! number of values in hm_RHM_RES(:) 
    ! HM input
    REAL(kind=dp) :: hm_rm
    REAL(kind=dp) :: hm_im
    REAL(kind=dp) :: hm_ll
    REAL(kind=dp) :: hm_rr
    INTEGER       :: hm_nangle
    REAL(kind=dp) :: hm_tilt
    ! HM output
    REAL(kind=dp) :: hm_RHM_RES(hm_nout)
    !                hm_rhm_res(1) = Number of scattered photons
    !                hm_rhm_res(2) = Incident energy flux
    !                hm_rhm_res(3) = scattered energy flux
    !                hm_rhm_res(4) = diffracted energy flux
    !                hm_rhm_res(5) = extinction cross section
    !                hm_rhm_res(6) = scattering cross section
    !                hm_rhm_res(7) = ssa
    !                hm_rhm_res(8) = asymetry factor
    REAL(kind=dp), allocatable :: hm_PGT(:,:,:)
    REAL(kind=dp), allocatable :: hm_thetaout(:)

    REAL(kind=dp), allocatable :: tmp_u(:)

    INTEGER :: i

    REAL(kind=dp) :: size_param
    REAL(kind=dp) :: r_eq
    REAL(kind=dp) :: v_crystal

    !========

    hm_ll = param_hm(1)
    hm_rr = param_hm(2)
    ! test ray tracing method validity
    ! compute the radius of the equivalent volume sphere
    v_crystal = (3.0_dp*SQRT(3.0_dp)/2.0_dp * hm_rr**2.0_dp * hm_ll)
    r_eq      = (3.0_dp/4.0_dp*v_crystal/xpi)**(1.0_dp/3.0_dp)
    size_param = 2.0_dp * xpi * r_eq / lamb
    if (size_param .lt. min_size_param_hm) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. CALL_RHM)  : ERROR '
       WRITE(*,*) '   For the ray tracing method to be valid, '
       WRITE(*,*) '   the size parameter must be greater than', min_size_param_hm
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   at lamb = ', lamb
       WRITE(*,*) '   2*PI*r_eq/lamb = ', size_param
       WRITE(*,*) '   r_eq= ', r_eq
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP
    endif

    hm_rm = REAL(refind)
    hm_im = AIMAG(refind)
    hm_tilt = param_hm(3)
    hm_nangle  = nang
    ALLOCATE(hm_PGT(4,4,hm_nangle),  &
         hm_thetaout(hm_nangle), &
         tmp_u(hm_nangle) )

    CALL hm_rhm1p1(nphot, lamb, hm_rm, hm_im,    &
         hm_ll, hm_rr, hm_tilt, hm_nangle,  & 
         hm_RHM_RES, hm_thetaout, hm_PGT)   

    if (hm_RHM_RES(7) .gt. 1.0_dp) THEN
       if (warning) then
          WRITE(*,*) ' '
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (in sub. call_RHM)     WARNING'
          WRITE(*,*) '   Scattering albedo can not be > 1.0 '
          WRITE(*,*) '   For particles ', ptype
          WRITE(*,*) '      w0 = ',hm_RHM_RES(7)
          WRITE(*,*) '   We reset it to 1.0'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' '
       endif
       hm_RHM_RES(7) = 1.0_dp
    endif

    ! sort in increasing u
    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(1,1,:))  ! F11
    call NM_SORT2(tmp_u, hm_PGT(1,1,:))  ! F11

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,2,:))  ! F22
    call NM_SORT2(tmp_u, hm_PGT(2,2,:))  ! F22

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,3,:))  ! F33
    call NM_SORT2(tmp_u, hm_PGT(3,3,:))  ! F33

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(4,4,:))  ! F44
    call NM_SORT2(tmp_u, hm_PGT(4,4,:))  ! F44

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,1,:))  ! F21
    call NM_SORT2(tmp_u, hm_PGT(2,1,:))  ! F21

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,4,:))  ! F34
    call NM_SORT2(tmp_u, hm_PGT(3,4,:))  ! F34


    ! write phase matrix in the common variable
    DO i= 1, hm_nangle
       u_phasemat(i)  = tmp_u(i)
       phasemat(1, i) = hm_PGT(1,1,i)               ! F11
       phasemat(2, i) = hm_PGT(2,2,i)*hm_PGT(1,1,i) ! F22
       phasemat(3, i) = hm_PGT(3,3,i)*hm_PGT(1,1,i) ! F33
       phasemat(4, i) = hm_PGT(4,4,i)*hm_PGT(1,1,i) ! F44
       phasemat(5, i) = hm_PGT(2,1,i)*hm_PGT(1,1,i) ! F21
       phasemat(6, i) = hm_PGT(3,4,i)*hm_PGT(1,1,i) ! F34
    ENDDO

    ! The computation is done over angle domains [thetamin:thetamax]
    ! so that we do not have 0 and 180 deg points in the final grid.
    ! We force the last and first angle to be
    ! 0 and 180 deg. 
    if(warning)then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_RHM)     WARNING'
       WRITE(*,*) ' The computation is done over angle domains [thetamin:thetamax]'
       WRITE(*,*) ' so that we do not have 0 and 180 deg points in the final grid.'
       WRITE(*,*) ' We force the last and first angle to be'
       WRITE(*,*) ' 0 and 180 deg. '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
    endif
    u_phasemat(1)         = -1.0_dp
    u_phasemat(hm_nangle) = 1.0_dp

    opt(1) = hm_rhm_res(5) ! Cext
    opt(2) = hm_rhm_res(7) ! ssa
    opt(3) = hm_rhm_res(8) ! asymetry factor

    DEALLOCATE(hm_PGT, hm_thetaout, tmp_u)

  END SUBROUTINE CALL_RHM

  !----------------------------------------------------------------

  SUBROUTINE CALL_IHM(nphot, param_hm, nang, refind, lamb, ptype, opt, phasemat, u_phasemat)

    USE MCONSTANTS, only : dp, &
         ptcle_nopt, &
         deg2rad, &
         hm_nmaxpar, &
         min_size_param_hm, &
         xpi, &
         max_len, &
         artdeco_path, &
         dir_libihm
    USE MCOMMON, ONLY : warning
    USE MUTILITY, only : trimcat
    !USE nr, ONLY : SORT2
    USE nm_nm, ONLY : NM_SORT2

    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nphot
    REAL(kind=dp), INTENT(IN)    :: param_hm(hm_nmaxpar)
    INTEGER, INTENT(IN)          :: nang
    COMPLEX(kind=dp), INTENT(IN) :: refind
    REAL(kind=dp), INTENT(IN)    :: lamb
    CHARACTER(LEN=max_len), INTENT(IN) :: ptype
    REAL(kind=dp), INTENT(OUT)   :: opt(ptcle_nopt)
    REAL(kind=dp), INTENT(OUT)   :: phasemat(6,nang)
    REAL(kind=dp), INTENT(OUT)   :: u_phasemat(nang)

    ! local variables
    ! HM input
    INTEGER, PARAMETER :: hm_nout = 8   ! number of values in hm_IHM_RES(:) 
    REAL(kind=dp) :: hm_rm
    REAL(kind=dp) :: hm_im
    REAL(kind=dp) :: hm_ll
    REAL(kind=dp) :: hm_rr
    INTEGER       :: hm_nangle
    REAL(kind=dp) :: hm_rm1
    REAL(kind=dp) :: hm_im1
    REAL(kind=dp) :: hm_rm2
    REAL(kind=dp) :: hm_im2
    REAL(kind=dp) :: hm_l
    REAL(kind=dp) :: hm_prop
    REAL(kind=dp) :: hm_reff_1
    REAL(kind=dp) :: hm_veff_1
    REAL(kind=dp) :: hm_reff_2
    REAL(kind=dp) :: hm_veff_2
    CHARACTER (LEN=max_len) :: hm_libdir

    ! HM output
    REAL(kind=dp) :: hm_IHM_RES(hm_nout)
    !                hm_ihm_res(1) = Number of scattered photons
    !                hm_ihm_res(2) = Incident energy flux
    !                hm_ihm_res(3) = scattered energy flux
    !                hm_ihm_res(4) = diffracted energy flux
    !                hm_ihm_res(5) = extinction cross section
    !                hm_ihm_res(6) = scattering cross section
    !                hm_ihm_res(7) = ssa
    !                hm_ihm_res(8) = asymetry factor
    REAL(kind=dp), allocatable :: hm_PGT(:,:,:)
    REAL(kind=dp), allocatable :: hm_thetaout(:)

    REAL(kind=dp), allocatable :: tmp_u(:)

    INTEGER :: i

    REAL(kind=dp) :: size_param
    REAL(kind=dp) :: r_eq
    REAL(kind=dp) :: v_crystal

    !========

    hm_ll = param_hm(1)
    hm_rr = param_hm(2)
    ! test ray tracing method validity
    ! compute the radius of the equivalent volume sphere
    v_crystal = (3.0_dp*SQRT(3.0_dp)/2.0_dp * hm_rr**2.0_dp * hm_ll)
    r_eq      = (3.0_dp/4.0_dp*v_crystal/xpi)**(1.0_dp/3.0_dp)
    size_param = 2.0_dp * xpi * r_eq / lamb
    if (size_param .lt. min_size_param_hm) then
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. CALL_IHM)  : ERROR '
       WRITE(*,*) '   For the ray tracing method to be valid, '
       WRITE(*,*) '   the size parameter must be greater than', min_size_param_hm
       WRITE(*,*) '   For particle : ', ptype
       WRITE(*,*) '   at lamb = ', lamb
       WRITE(*,*) '   2*PI*r_eq/lamb = ', size_param
       WRITE(*,*) '   r_eq= ', r_eq
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP
    endif

    hm_rm     = REAL(refind)  ! real part of  particles refractive index
    hm_im     = AIMAG(refind) ! imaginary part of particle refractive index 

    ! INCLUSIONS 1 are taken to be air bubbles :
    ! the real part of air refractive index is taken to be 1.0 at any wavelength
    ! so that the real part of the air inclusions refractive index surrounded
    ! by ice will be 1/n_ice
    ! We assume the air inclusions to be non absorbing at all wavelengths k = 0.0
    ! the refractive index of air inclusions regarding surronding medium is then 1 / m_ice
    hm_rm1    = 1.0_dp / REAL(refind)   ! real part of relative refractive index for inclusion type 1 
    hm_im1    = 0.0_dp                  ! imaginary part of refractive index for inclusion type 1
    hm_rm2    = 1.00_dp      ! real part of relative refractive index for inclusion type 2 
    hm_im2    = 0.00_dp      ! imaginary part refractive index for inclusion type 2 
    hm_l      = param_hm(3) ! mean free path between inclusion
    hm_prop   = 0.0_dp                     ! N_inclusion_2 / N_inclusion_1
    hm_reff_1 = param_hm(4)
    hm_veff_1 = param_hm(5)
    hm_reff_2 = hm_reff_1
    hm_veff_2 = hm_veff_1
    hm_nangle = nang
    hm_libdir = TRIMCAT(artdeco_path, dir_libihm)

    ALLOCATE(hm_PGT(4,4,hm_nangle),  &
         hm_thetaout(hm_nangle), &
         tmp_u(hm_nangle) )

    CALL hm_ihm2p0(nphot,&
         lamb, &
         hm_rm, &
         hm_im, &
         hm_ll, &
         hm_rr, &
         hm_rm1, &
         hm_im1, &
         hm_rm2, &
         hm_im2, &
         hm_l, &
         hm_prop, &
         hm_reff_1, &
         hm_veff_1, &
         hm_reff_2, &
         hm_veff_2, &
         hm_nangle, &
         hm_libdir, &
         hm_IHM_RES, &
         hm_thetaout, &
         hm_PGT)   

    if (hm_IHM_RES(7) .gt. 1.0_dp) THEN
       if(warning) then
          WRITE(*,*) ' '
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (in sub. call_IHM)     WARNING'
          WRITE(*,*) '   Scattering albedo can not be > 1.0 '
          WRITE(*,*) '   For particles ', ptype
          WRITE(*,*) '      w0 = ',hm_IHM_RES(7)
          WRITE(*,*) '   We reset it to 1.0'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' '
       endif
       hm_IHM_RES(7) = 1.0_dp
    endif

    ! sort in increasing u
    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(1,1,:))  ! F11
    call NM_SORT2(tmp_u, hm_PGT(1,1,:))  ! F11

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,2,:))  ! F22
    call NM_SORT2(tmp_u, hm_PGT(2,2,:))  ! F22

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,3,:))  ! F33
    call NM_SORT2(tmp_u, hm_PGT(3,3,:))  ! F33

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(4,4,:))  ! F44
    call NM_SORT2(tmp_u, hm_PGT(4,4,:))  ! F44

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(2,1,:))  ! F21
    call NM_SORT2(tmp_u, hm_PGT(2,1,:))  ! F21

    tmp_u(:) = cos(hm_thetaout(:)*deg2rad)
    !call sort2(tmp_u, hm_PGT(3,4,:))  ! F34
    call NM_SORT2(tmp_u, hm_PGT(3,4,:))  ! F34


    ! write phase matrix in the common variable
    DO i= 1, hm_nangle
       u_phasemat(i)  = tmp_u(i)
       phasemat(1, i) = hm_PGT(1,1,i)               ! F11
       phasemat(2, i) = hm_PGT(2,2,i)*hm_PGT(1,1,i) ! F22
       phasemat(3, i) = hm_PGT(3,3,i)*hm_PGT(1,1,i) ! F33
       phasemat(4, i) = hm_PGT(4,4,i)*hm_PGT(1,1,i) ! F44
       phasemat(5, i) = hm_PGT(2,1,i)*hm_PGT(1,1,i) ! F21
       phasemat(6, i) = hm_PGT(3,4,i)*hm_PGT(1,1,i) ! F34
    ENDDO

    ! The computation is done over angle domains [thetamin:thetamax]
    ! so that we do not have 0 and 180 deg points in the final grid.
    ! We force the last and first angle to be
    ! 0 and 180 deg.
    if (warning) then 
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_IHM)     WARNING'
       WRITE(*,*) ' The computation is done over angle domains [thetamin:thetamax]'
       WRITE(*,*) ' so that we do not have 0 and 180 deg points in the final grid.'
       WRITE(*,*) ' We force the last and first angle to be'
       WRITE(*,*) ' 0 and 180 deg. '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
    endif
    u_phasemat(1)         = -1.0_dp
    u_phasemat(hm_nangle) =  1.0_dp

    opt(1) = hm_ihm_res(5) ! Cext
    opt(2) = hm_ihm_res(7) ! ssa
    opt(3) = hm_ihm_res(8) ! asymetry factor

    DEALLOCATE(hm_PGT, hm_thetaout, tmp_u)

  END SUBROUTINE CALL_IHM

END MODULE MGET_OPT
