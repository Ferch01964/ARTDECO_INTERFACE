

MODULE MF2PY_UTILITY

  USE MCONSTANTS, only : max_len

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: print_kdis_species, &
       set_kdis_species, &
       allocate_ptcle_char_arr, &
       set_ptcle_vdist, &
       alloc_ptcle_betal, &
       print_ptcle_betal, &
       allocate_gas_u_atm, &
       set_gas_u_atm

CONTAINS 

  subroutine print_ptcle_betal(i, j, k, l)

    USE MCOMMON, only : ptcle_betal
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: k
    INTEGER, INTENT(IN) :: l

    print*, "ptcle_betal = ",  ptcle_betal(i, j, k, l)
    
  end subroutine print_ptcle_betal

  !----------------------------------------------------------------
 
 subroutine alloc_ptcle_betal
    ! because I do not know how to make a dimension
    ! starting at 0 from ptyhon side allocation

    USE MCOMMON, only : ptcle_betal, nptcle, nmaxbetal, nlambda
    USE MCONSTANTS, only : undef_dp

    IMPLICIT NONE

    ALLOCATE(ptcle_betal(nptcle, 6, 0:nmaxbetal, nlambda))
    ptcle_betal(:, :, :, :) = undef_dp
    ptcle_betal(:, :, 0, :) = 0.0
    
  end subroutine alloc_ptcle_betal

  !----------------------------------------------------------------

  subroutine set_ptcle_vdist(indarr, vdist)

    USE MCOMMON, only : ptcle_vdist_type

    IMPLICIT NONE
    CHARACTER (LEN=max_len), INTENT(IN) :: vdist
    INTEGER, INTENT(IN) :: indarr
  
    ptcle_vdist_type(indarr) = vdist
    
  end subroutine set_ptcle_vdist

  !----------------------------------------------------------------

  subroutine allocate_ptcle_char_arr(nptcle)

    USE MCONSTANTS, only : undef_c
    USE MCOMMON, only : ptcle_type,ptcle_vdist_type

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nptcle

    ALLOCATE(ptcle_type(nptcle), ptcle_vdist_type(nptcle))
    
    ptcle_type(:)                 = undef_c
    ptcle_vdist_type(:)           = undef_c
    
  end subroutine allocate_ptcle_char_arr

  !----------------------------------------------------------------

  subroutine allocate_gas_u_atm(n)

    USE MCONSTANTS, only : undef_c
    USE MCOMMON, only : gas_u_atm

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    ALLOCATE(gas_u_atm(n))
    gas_u_atm(:)                 = undef_c
    
  end subroutine allocate_gas_u_atm

  !----------------------------------------------------------------

  subroutine set_gas_u_atm(indarr, name)

    USE MCOMMON, only : gas_u_atm

    IMPLICIT NONE
    CHARACTER (LEN=max_len), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: indarr
  
    gas_u_atm(indarr) = name
    
  end subroutine set_gas_u_atm

  !----------------------------------------------------------------

  function print_kdis_species(indarr)

    USE MCOMMON, only: kdis_nsp, kdis_species

    IMPLICIT NONE
    CHARACTER (LEN=max_len) :: print_kdis_species
    INTEGER, INTENT(IN) :: indarr

    print_kdis_species = kdis_species(indarr)

  end function print_kdis_species

  !----------------------------------------------------------------

  subroutine set_kdis_species(indarr, spec)

    USE MCOMMON, only: kdis_nsp, kdis_species

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: indarr
    CHARACTER (LEN=max_len) :: spec

    kdis_species(indarr) = spec

  end subroutine set_kdis_species


END MODULE MF2PY_UTILITY
