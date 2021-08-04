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

  SUBROUTINE compute_dummy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    DO dum_int = 1,10

      dummy_ar( dum_int ) = 1.0D0

    END DO

  END

  SUBROUTINE compute_shell_grid_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get data on the shell grid - list of |\K
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    DO ind_1 = 1, gr1_size
      shell_en_XS1( ind_1 ) = CDABS( v_x(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two + &
                              CDABS( v_y(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two + &
                              CDABS( v_z(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two
      shell_es_XS1( ind_1 ) = CDABS( w_ux(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two + &
                              CDABS( w_uy(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two + &
                              CDABS( w_uz(sh_gr1_x,sh_gr1_y,sh_gr1_z ) )**two
    END DO
    DO ind_2 = 1, gr2_size
      shell_en_XS2( ind_2 ) = CDABS( v_x(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two + &
                              CDABS( v_y(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two + &
                              CDABS( v_z(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two
      shell_es_XS2( ind_2 ) = CDABS( w_ux(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two + &
                              CDABS( w_uy(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two + &
                              CDABS( w_uz(sh_gr2_x,sh_gr2_y,sh_gr2_z ) )**two
    END DO

    CALL write_shell_grid_data
    ! REF-> <<< system_advoutput >>>

  END

END MODULE system_advfunctions
