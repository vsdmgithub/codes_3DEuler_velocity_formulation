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

  SUBROUTINE find_k_thermalising
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the k_th(t)
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::k_ind

    en_min = spectral_energy_avg( k_th )
    k_min  = k_th

    DO k_ind = 2, k_G - 1

      IF( spectral_energy_avg( k_ind )  .LT. en_min  ) THEN

        k_min  = k_ind
        en_min = spectral_energy_avg( k_min )

      END IF

    END DO

    k_th = k_min

  END

  SUBROUTINE find_spectral_thermalising_coefficient
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get a coefficient of thermalisation in log scale for last few spectra
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::k_ind

    th_supp_factor = one

    DO k_ind = k_iner_ref + 1 , k_th

      en_ref   = spectral_energy_avg( k_iner_ref ) * ( DBLE( k_ind / k_iner_ref ) ** ( -fivthird ) )

      IF ( spectral_energy_avg( k_ind ) .GT. tol ) THEN

        th_coeff = DLOG( spectral_energy_avg( k_ind ) / en_ref )

        IF ( th_coeff .GT. th_coeff_threshold ) THEN

          th_supp_factor( k_ind ) = DEXP( -th_coeff / two )

        END IF

      END IF

    END DO

  END

  SUBROUTINE thermal_suppression_filter
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the velocity suppressed from thermalisation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::k_ind

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

        k_ind = shell_no( i_x, i_y, i_z )

  			IF ( th_supp_factor( k_ind ) .LT. one ) THEN

          v_x( i_x, i_y, i_z ) = v_x ( i_x, i_y, i_z ) * th_supp_factor( k_ind )
          v_y( i_x, i_y, i_z ) = v_y ( i_x, i_y, i_z ) * th_supp_factor( k_ind )
          v_z( i_x, i_y, i_z ) = v_z ( i_x, i_y, i_z ) * th_supp_factor( k_ind )

  			END IF

      END IF

		END DO
		END DO
		END DO

  END

  SUBROUTINE thermal_cutoff_filter
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the spectrum with cutoff at k_th
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_th_2 = DBLE(k_th * k_th)

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

			IF ( k_2( i_x, i_y, i_z ) .GT. k_th_2 ) THEN

        v_x( i_x, i_y, i_z ) = c0
        v_y( i_x, i_y, i_z ) = c0
        v_z( i_x, i_y, i_z ) = c0

			END IF

		END DO
		END DO
		END DO

  END

END MODULE system_advfunctions
