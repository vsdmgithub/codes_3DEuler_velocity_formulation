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
  USE system_VTR
  USE system_VTK

  IMPLICIT NONE

  TYPE(VTR_file_handle)::fd
  ! creating dataypes to store paraview files for 3d viewing

  CONTAINS

  SUBROUTINE write_PVD_shell_grid
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in PVD format
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    ! WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) &
                // 'modal_t'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  SHELL DATA
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=file_name,FD=fd)

    CALL  VTR_write_mesh(FD=fd,X=k_x_axis,Y=k_y_axis,Z=k_z_axis)

    shell_en_matrix                                   = zero
    shell_es_matrix                                   = zero

    ! DO ind_1 = 1, gr1_size
    !
    !   IF ( shell_en_XS1( ind_1 ) .GT. tol ) THEN
    !     shell_en_matrix( sh_gr1_x( ind_1 ), sh_gr1_y( ind_1 ), sh_gr1_z( ind_1 ) ) = DLOG( shell_en_XS1( ind_1 ) / tol )
    !   END IF
    !   IF ( shell_es_XS1( ind_1 ) .GT. tol ) THEN
    !     shell_es_matrix( sh_gr1_x( ind_1 ), sh_gr1_y( ind_1 ), sh_gr1_z( ind_1 ) ) = DLOG( shell_es_XS1( ind_1 ) / tol )
    !   END IF
    !
    ! END DO
    ! !
    ! DO ind_2 = 1, gr2_size
    !
    !   IF ( shell_en_XS2( ind_2 ) .GT. tol ) THEN
    !     shell_en_matrix( sh_gr2_x( ind_2 ), sh_gr2_y( ind_2 ), sh_gr2_z( ind_2 ) ) = DLOG( shell_en_XS2( ind_2 ) / tol )
    !   END IF
    !   IF ( shell_es_XS2( ind_2 ) .GT. tol ) THEN
    !     shell_es_matrix( sh_gr2_x( ind_2 ), sh_gr2_y( ind_2 ), sh_gr2_z( ind_2 ) ) = DLOG( shell_es_XS2( ind_2 ) / tol )
    !   END IF
    !
    ! END DO

    DO i_x = kMin_x, kMax_x
    DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z
      energy_mode = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                    CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                    CDABS( v_z( i_x, i_y, i_z ) ) ** two
      IF ( hf * energy_mode .GT. tol ) THEN
        shell_en_matrix( i_x, i_y, i_z ) = DLOG( hf * energy_mode / tol )
      END IF
      enstrophy_mode = CDABS( w_vx( i_x, i_y, i_z ) ) ** two + &
                       CDABS( w_vy( i_x, i_y, i_z ) ) ** two + &
                       CDABS( w_vz( i_x, i_y, i_z ) ) ** two
      IF ( hf * enstrophy_mode .GT. tol ) THEN
        shell_es_matrix( i_x, i_y, i_z ) = DLOG( hf * enstrophy_mode / tol )
      END IF

    END DO
    END DO
    END DO

    CALL VTR_write_var(FD=fd, NAME='Energy',FIELD=shell_en_matrix)

    CALL VTR_write_var(FD=fd, NAME='Enstrophy',FIELD=shell_es_matrix)

    CALL  VTR_close_file(FD=fd)

    ! CALL  VTR_collect_file( FD = fd )
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_shell_grid_XS
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the location of shell grid surface
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'sh_gr1_XS'// '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 4230, FILE = file_name )
    DO ind_1 = 1, gr1_size
      !
      WRITE(4230,f_d8p4,ADVANCE ='no')  k_x_axis( sh_gr1_x( ind_1 ) )
      WRITE(4230,f_d8p4,ADVANCE ='no')  k_y_axis( sh_gr1_y( ind_1 ) )
      WRITE(4230,f_d8p4,ADVANCE ='yes') k_z_axis( sh_gr1_z( ind_1 ) )
      !
    END DO

    CLOSE(4230)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'sh_gr2_XS'// '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 4232, FILE = file_name )
    DO ind_2 = 1, gr2_size

      WRITE(4232,f_d8p4,ADVANCE ='no')  k_x_axis( sh_gr2_x( ind_2 ) )
      WRITE(4232,f_d8p4,ADVANCE ='no')  k_y_axis( sh_gr2_y( ind_2 ) )
      WRITE(4232,f_d8p4,ADVANCE ='yes') k_z_axis( sh_gr2_z( ind_2 ) )

    END DO

    CLOSE(4232)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_shell_grid_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the data of points on the shell  of given k
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'sh_gr1_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 4231, FILE = file_name )
    DO ind_1 = 1, gr1_size

      WRITE(4231,f_d32p17,ADVANCE ='no')   shell_en_XS1( ind_1 )
      WRITE(4231,f_d32p17,ADVANCE ='yes')  shell_es_XS1( ind_1 )

    END DO

    CLOSE(4231)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'sh_gr2_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 4233, FILE = file_name )
    DO ind_2 = 1, gr2_size

      WRITE(4233,f_d32p17,ADVANCE ='no')   shell_en_XS2( ind_2 )
      WRITE(4233,f_d32p17,ADVANCE ='yes')  shell_es_XS2( ind_2 )

    END DO

    CLOSE(4233)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
