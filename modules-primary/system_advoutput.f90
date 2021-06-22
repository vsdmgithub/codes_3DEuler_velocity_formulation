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

  SUBROUTINE write_vx_dot_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_dot_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 456, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_y = Nh - 1
    DO i_x = 0, Nh
    DO i_z = 0, Nh
         WRITE(456,f_d32p17,ADVANCE='no')  str_xy(i_x, i_y, i_z) + w_uz( i_x, i_y, i_z)
         WRITE(456,f_d32p17,ADVANCE='no')  two * str_yy(i_x, i_y, i_z)
         WRITE(456,f_d32p17,ADVANCE='yes') str_yz(i_x, i_y, i_z) + w_ux(i_x, i_y, i_z)
    END DO
    END DO

    CLOSE(456)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vx_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 455, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_y = Nh - 1
    DO i_x = 0, Nh
    DO i_z = 0, Nh
         WRITE(455,f_d32p17,ADVANCE='no')  w_ux(i_x, i_y, i_z)
         WRITE(455,f_d32p17,ADVANCE='no')  w_uy(i_x, i_y, i_z)
         WRITE(455,f_d32p17,ADVANCE='yes') w_uz(i_x, i_y, i_z)
    END DO
    END DO

    CLOSE(455)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vx_stretching_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_alp_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 788, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! vx_stretching_max = MAXVAL( DABS( vx_stretching ) )
    ! For Normalization by the infinite norm of the vx stretching factor.

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         WRITE(788,f_d32p17,ADVANCE='yes') vx_stretching(i_x,i_y,i_z)
    END DO
    END DO

    CLOSE(788)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_loc_vx_stretching_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_dot_loc_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 789, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! loc_vx_stretching_max = MAXVAL( DABS( vx_stretching - bck_vx_stretching ) )
    ! For Normalization by the infinite norm of the vx stretching factor.

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         loc_stretching = vx_stretching(i_x,i_y,i_z) - bck_vx_stretching(i_x,i_y,i_z)
         WRITE(789,f_d32p17,ADVANCE='no')  vx_stretching( i_x, i_y, i_z )
         WRITE(789,f_d32p17,ADVANCE='yes') loc_stretching
    END DO
    END DO

    CLOSE(789)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
