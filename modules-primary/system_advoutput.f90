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

  SUBROUTINE write_k_th_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   K _ T H     V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'k_th.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4304, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(4304,f_d8p4,ADVANCE   ='no')  time_now
    WRITE(4304,f_i6,ADVANCE ='yes') k_th

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4304)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_thermalising_coefficient
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the spectral thermalising coefficient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) &
                // 'spectral_coeff_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1201, FILE = file_name )
    DO k_no = k_iner_ref , k_G

      WRITE(1201,f_i8,ADVANCE  ='no')       k_no
      WRITE(1201,f_d32p17,ADVANCE ='yes')   th_supp_factor( k_no )

    END DO
    CLOSE(1201)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END
  
END MODULE system_advoutput
