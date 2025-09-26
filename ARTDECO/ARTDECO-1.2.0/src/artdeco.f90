

PROGRAM ARDECO

  USE MCONSTANTS, only : max_len
  USE MCOMMON, only : rt_model, &
       nptcle, &
       opt_only, &
       do_rt, &
       mode, &
       artdeco_infile_name

  USE MINOUT
  USE MGET_OPT
  USE MGET_BETAL
  USE MSOLRAD
  USE MLAYERS
  USE MCALL_AD
  USE MCALL_OD
  USE MCALL_DOAD
  USE MCALL_MCRAD1D
  USE MSINSCA
  USE MBAND_INT

  IMPLICIT NONE

  character(len=max_len) :: progname

  !===================================

  ! Get program name
  call get_command_argument(0, progname)
  progname=(progname(index(progname,'/',.true.)+1:))

  if (command_argument_count() < 1) then
    write(0,*) 'Syntax  : ',trim(progname),'   infile_name'
    stop
  endif

  call get_command_argument(1,artdeco_infile_name)
  
  ! must be the first called subroutine
  CALL READ_DATA

  IF (nptcle .gt. 0) THEN
     CALL GET_OPT
     IF (opt_only .eqv. .FALSE.) CALL GET_BETAL
  ENDIF

  IF (do_rt) THEN

     CALL GET_SOLRAD

     CALL GET_LAYERS

     SELECT CASE (rt_model)
     CASE('addoub')
        CALL CALL_AD
     CASE('doad')
        CALL CALL_DOAD
     CASE('disort')
        CALL CALL_OD
     CASE('mcrad1d')
        CALL CALL_MCRAD1D
     CASE('sinsca')
        CALL CALL_SINSCA
     END SELECT

     if (mode.ne.'mono') CALL GET_BAND_INT

  ENDIF

  ! must be the last called subroutine
  CALL WRITE_DATA

END PROGRAM ARDECO


