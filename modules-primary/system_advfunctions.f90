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

  SUBROUTINE compute_shell_grid_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get data on the shell grid - list of |\K
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN

      CALL write_shell_grid_XS
      ! REF-> <<< system_advoutput >>>

    END IF
    
    DO ind_1 = 1, gr1_size
      shell_en_XS1( ind_1 ) = CDABS( v_x(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two + &
                              CDABS( v_y(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two + &
                              CDABS( v_z(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two
      shell_es_XS1( ind_1 ) = CDABS( w_vx(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two + &
                              CDABS( w_vy(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two + &
                              CDABS( w_vz(sh_gr1_x(ind_1),sh_gr1_y(ind_1),sh_gr1_z(ind_1) ) )**two
    END DO
    DO ind_2 = 1, gr2_size
      shell_en_XS2( ind_2 ) = CDABS( v_x(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two + &
                              CDABS( v_y(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two + &
                              CDABS( v_z(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two
      shell_es_XS2( ind_2 ) = CDABS( w_vx(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two + &
                              CDABS( w_vy(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two + &
                              CDABS( w_vz(sh_gr2_x(ind_2),sh_gr2_y(ind_2),sh_gr2_z(ind_2) ) )**two
    END DO

    CALL write_shell_grid_data
    ! REF-> <<< system_advoutput >>>

  END

END MODULE system_advfunctions
