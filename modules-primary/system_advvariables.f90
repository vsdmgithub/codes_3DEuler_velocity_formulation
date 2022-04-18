! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------

! ##################
! MODULE: system_advvariables
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED VARIABLES TO DO ANALYSIS IN 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_basicfunctions

  IMPLICIT NONE
  ! _________________________
  ! VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND =4) ::k_min
  INTEGER(KIND =4) ::k_th
  INTEGER(KIND =4) ::k_iner_ref
  INTEGER(KIND =4) ::t_step_supp
  DOUBLE PRECISION ::en_min,en_ref
  DOUBLE PRECISION ::k_th_2
  DOUBLE PRECISION ::time_supp
  DOUBLE PRECISION ::th_coeff,th_coeff_threshold
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::th_supp_factor

  CONTAINS

  SUBROUTINE allocate_thermal_filters
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( th_supp_factor(0:max_shell_no) )

    k_iner_ref     = FLOOR( DBLE( N ) / fiv )
    ! Reference for inertial range wavenumber

    th_coeff_threshold = 1.1D0
    ! ratio of spectral energy to reference energy with fluctuations included

    k_th           = k_G - 1
    ! Thermalisation wavenumber initial status

    time_supp      = ( DBLE( 2.0D0 / k_th ) ** twothird ) * ( length / ( six * v_rms_1D ) ) * 0.75D0
    ! Time scale for the thermalising wavenumber

    CALL time_to_step_convert(time_supp,t_step_supp,dt)
    ! REF-> <<< system_auxilaries >>>

  END

END MODULE system_advvariables
