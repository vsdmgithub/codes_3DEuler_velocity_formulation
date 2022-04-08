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
! MODULE: system_advfunctions
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED FUNCTIONS MODULE FOR 3D EULER ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advfunctions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicfunctions
  USE system_advvariables
  USE system_advoutput

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE find_k_thermal
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the k_th(t)
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::k_min,k_ind
    DOUBLE PRECISION::en_min

    en_min = spectral_energy_avg( k_G - 1 )
    k_min  = k_G - 1

    DO k_ind = 5, k_G - 1

      IF( spectral_energy_avg( k_ind )  .LT. en_min  ) THEN

        k_min = k_ind

      END IF

    END DO

    k_th = k_min

  END

  SUBROUTINE thermal_suppression_filter
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the thermally suppressed velocity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_th_2 = DBLE(k_th * k_th)

	    DO i_x = 0, Nh
	    DO i_y = -Nh, Nh - 1
	    DO i_z = -Nh, Nh - 1

				IF ((k_2(i_x,i_y,i_z)) .GT. k_th_2 ) THEN

          v_x( i_x, i_y, i_z ) = c0
          v_y( i_x, i_y, i_z ) = c0
          v_z( i_x, i_y, i_z ) = c0

				END IF

			END DO
			END DO
			END DO

  END

END MODULE system_advfunctions
