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
! PROGRAM THAT SOLVES 3D INCOMPRESSIBLE EULER EQUATION
! LAST MODIFIED: 3 JUNE 2021
! ##################

! ##################
! LIST OF MODULES:
! ___________________
! PRIMARY
! ===================
! system_main
! system_advectionsolver
! system_vorticitysolver
! system_basicfunctions
! system_advfunctions
! system_initialcondition
! system_basicvariables
! system_advvariables
! system_constants
! system_basicoutput
! system_advoutput
! system_auxilaries
! system_fftw
! system_pvd_output
! system_timer
! -------------------
! ##################

PROGRAM euler
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the work is done in the modules. Calling a few would finish
! the code.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
!  MODULES
!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
USE system_main
USE system_test
USE system_timer

	IMPLICIT NONE
	! _________________________
	! LOCAL VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!

  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  I  N  I  T  I  A  L  I  Z  A  T  I  O  N
  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL read_input
	! REF-> <<< system_variables >>>

  CALL init_global_variables
	! REF-> <<< system_variables >>>

	CALL init_global_arrays
	! REF-> <<< system_variables >>>

	! We get all the variables, and global arrays ready allocated

	test_code		= 'n'
	! Simple way to ON or OFF the testing of the simulation - Measures the time for evolution

  run_code 		= 'y'
	! Simple way to ON or OFF the running of the simulation.

  IF ( run_code .EQ. 'y' ) THEN

		CALL start_run_timer
		! Clocks the date and time - PROGRAM STARTS
		! REF-> <<< system_timer >>>

    CALL pre_analysis
    ! Allocating the evolution arrays, if everything is set, 'check_status' will be 1.
		! REF-> <<< system_main >>>

    IF ( check_status .EQ. 1 ) THEN

			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('   S  I  M  U  L  A  T  I  O  N        S  T  A  R  T  S '))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,*)

      CALL time_evolution
      ! Solve the 3D Euler equation, in discrete time using pseudospectral method.
			! REF-> <<< system_main >>>

			IF ( debug_error .NE. 1 ) THEN

			  CALL post_analysis
	      ! Does the post-analysis, final outputs and deallocation
				! REF-> <<< system_main >>>

			END IF

    END IF

    IF ( state_sim .EQ. 1 ) THEN

			WRITE(*,*)
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('    S  I  M  U  L  A  T  I  O  N        E  N  D  S  '))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))

    END IF

		CALL end_run_timer
		! Clocks the date and time - PROGRAM STARTS
		! REF-> <<< system_timer >>>

  END IF

	IF ( test_code .EQ. 'y' ) THEN

		CALL pre_analysis
    ! Allocating the evolution arrays
		! REF-> <<< system_main >>>

		CALL test_fft_time
		! Check forward and inverse fft time
		! REF-> <<< system_test >>>

		CALL test_evolution_time
		! Checks time for one evolution step - predicts total time
		! REF-> <<< system_test >>>

	END IF

END program euler
