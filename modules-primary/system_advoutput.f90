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
! MODULE: system_advoutput
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED OUTPUT MODULE - RELATED TO ADVANCED FUNCTION MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all output calls from advanced functions module.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advvariables
  USE system_basicoutput

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE write_dummy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the list of exponents
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! file_name = TRIM( ADJUSTL( file_address ) ) // 'dummy.dat'
    !
    ! !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! OPEN(unit = 8880, file = file_name )
    !
    ! DO dum_int = 1,10
    !
    !   WRITE(8880,f_d8p4,ADVANCE ='yes') dummy_ar( dum_int )
    !
    ! END DO
    !
    ! !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! CLOSE(8880)

  END

  SUBROUTINE write_adv_simulation_details
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   Write the details of the simulation - only the advanced parameter
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    file_name = TRIM(ADJUSTL(file_address))//'system_details_adv'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  I  M  U  L  A  T  I  O  N     D  E  T  A  I  L  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(UNIT =433,FILE=TRIM( ADJUSTL( file_name ) ) // '.dat')

    WRITE(433,"(A50)")TRIM(ADJUSTL('-------------------------------------------------------'))
    WRITE(433,"(A50)")TRIM(ADJUSTL('------3D EULER EQUATION (INCOMPRESSIBLE)--------------'))
    WRITE(433,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))
    WRITE(433,*)
    WRITE(433,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(433,"(A50)")TRIM(ADJUSTL('---------ADV PARAMETERS OF SIMULATION----------'))
    WRITE(433,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(433,*)
    WRITE(433,"(A30,A2,ES8.2)")  'Thermal filter timestep','= ',time_supp
    WRITE(433,"(A30,A2,I8)")     'Filtering interval ','= ',t_step_supp
    WRITE(433,"(A30,A2,ES8.2)")  'Thermal filtering coeff','= ',th_coeff_threshold
    WRITE(433,"(A30,A2,I8)")     'k_inertial ref ','= ',k_iner_ref
    WRITE(433,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))

    CLOSE(433)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
