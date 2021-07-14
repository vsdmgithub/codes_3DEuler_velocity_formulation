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
! MODULE: system_basicoutput
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC OUTPUT MODULE - RELATED TO BASIC FUNCTION MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all basic outputs produced in the simulation.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables

  IMPLICIT NONE
  ! _________________________
  ! OUTPUT VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(LEN =140)::file_name
  CHARACTER(LEN =40) ::file_time
  CHARACTER(LEN =40) ::path_dir
  CHARACTER(LEN =80) ::type_sim
  CHARACTER(LEN =40) ::name_sim
  CHARACTER(LEN =100)::file_address
  CHARACTER(LEN =40) ::sub_dir_3D
  CHARACTER(LEN =40) ::sub_dir_2D
  CHARACTER(LEN =40) ::sub_dir_sp
  CHARACTER(LEN =40) ::sub_dir

  CONTAINS

  SUBROUTINE prepare_output
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Name the folders, create them, open files to write system_basicoutput .
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    CALL name_output_dir
    ! Names all the directories where output is stored

    CALL create_output_dir
    ! Creates the directories

    CALL write_simulation_details
    ! Writes the parameters used in the simulation

	END

  SUBROUTINE name_output_dir
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Name the folders, file address
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    path_dir    =   '../euler_data/'
    ! path of the main directory relative to this file.

    sub_dir_3D  =   '3D_data/'
    ! Sub directory name to store 3D data - large file sizes.

    sub_dir_2D  =   '2D_data/'
    ! Sub directory name to store section files (2D data)

    sub_dir_sp  =   'k_data/'
    ! Sub directory name to store spectral data

    ! type_sim    =   TRIM( ADJUSTL( N_char ) ) // '/'
    type_sim    =   'VX_S_' // TRIM( ADJUSTL( N_char ) ) // '/'
    ! type of simulation, the data is storing

    CALL get_simulation_name(name_sim)
    ! Creating dated and timed name for the simulation for this particular type

    ! name_sim    =   'test_sim'
    ! Use this to give CUSTOM SIMULATION NAME

    file_address =   TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) //  &
                     TRIM( ADJUSTL( name_sim ) ) // '/'
    ! Address should be added to all file names, if needed sub-dir can be declared later and appended to
    ! this address

	END

  SUBROUTINE create_output_dir
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Create the directories.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) )

    CALL SYSTEM('mkdir '// TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) )

    ! Command to create the main directory and sub directories (name_sim) in the desired path
    ! If exists already, it won't be an error

	END

  SUBROUTINE write_test_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the test_data of any new verification/validation conducted in the code
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::test_file_name

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      test_file_name = TRIM( ADJUSTL( file_address ) ) // 'test_data.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 9009, file = test_file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(9009,f_d8p4,ADVANCE   ='no')   time_now
    WRITE(9009,f_d32p17,ADVANCE ='yes')  u_x( 2, 2, 2)

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(9009)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_simulation_details
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   Write the details of the simulation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    file_name = TRIM(ADJUSTL(file_address))//'system_details'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  I  M  U  L  A  T  I  O  N     D  E  T  A  I  L  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(UNIT =233,FILE=TRIM( ADJUSTL( file_name ) ) // '.dat')

    WRITE(233,"(A50)")TRIM(ADJUSTL('-------------------------------------------------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('------3D EULER EQUATION (INCOMPRESSIBLE)--------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(233,"(A20,A2,I3,A2,I3,A2,I3)") 'L x W x H ','= ',N_x,'X',N_y,'X',N_z
    WRITE(233,"(A20,A2,I8)") 'Trunc. Mode ',          '= ',k_G
    WRITE(233,"(A20,A2,ES8.2)") 'Time step ',         '= ',dt
    WRITE(233,"(A20,A2,I8)") 'CFL ratio ',            '= ',CFL_system
    WRITE(233,"(A20,A2,I8)") 'Total time steps ',     '= ',t_step_total
    WRITE(233,"(A20,A2,F8.4)") 'Total time ',         '= ',time_total
    WRITE(233,"(A20,A2,I8)") 'No of saves ',          '= ',no_of_saves
    WRITE(233,"(A20,A2,I8)") 'No of PVD saves ',      '= ',no_of_PVD_saves
    WRITE(233,"(A20,A2,F8.4)") 'Initial energy ',     '= ',energy
    WRITE(233,"(A20,A2,F8.4)") 'Initial enstrophy ',  '= ',enstrophy
    WRITE(233,"(A20,A2,ES8.2)") 'Initial compress ',  '= ',k_dot_v_norm
    WRITE(233,"(A20,A2,A8)") 'Initial cond',          '= ',TRIM( ADJUSTL( IC_type ) )
    WRITE(233,"(A20,A2,I8)") 'Total modes ',          '= ',tot_modes
    WRITE(233,"(A20,A2,I8)") 'Active modes',          '= ',tot_active_modes
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))

    CLOSE(233)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) &
                // 'no_modes_in_shell.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  M  O  D  E  S    I  N     S  H  E  L  L
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1991, FILE = file_name )
    DO k_no = 1 , max_shell_no

      WRITE(1991,f_i4,ADVANCE  ='no')                    k_no
      WRITE(1991,f_i8,ADVANCE ='yes') count_modes_shell( k_no )

    END DO
    CLOSE(1991)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END

  SUBROUTINE write_temporal_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'energy_vs_time.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4004, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(4004,f_d8p4,ADVANCE   ='no')  time_now
    WRITE(4004,f_d32p17,ADVANCE ='no')  energy
    WRITE(4004,f_d32p17,ADVANCE ='yes') enstrophy

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4004)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the spectral energy into a file with time index
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) &
                // 'spectral_energy_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1001, FILE = file_name )
    DO k_no = 1 , max_wave_no

      WRITE(1001,f_i8,ADVANCE  ='no')                            k_no
      WRITE(1001,f_d32p17,ADVANCE ='no')    spectral_energy(     k_no )
      WRITE(1001,f_d32p17,ADVANCE ='yes')   spectral_energy_avg( k_no )

    END DO
    CLOSE(1001)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_spectral_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a spectral velocity, in such a way that it can be used
  ! for input in the next simulation using 'IC_from_file' subroutine.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_3D )) &
                // 'spectral_velocity_' //TRIM( ADJUSTL( N_char ) ) // '_input.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where energy vs time will be written. With additional data

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( unit = 73, file = file_name )

    DO i_x = kMin_x, kMax_x
    DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

      WRITE(73,f_c32p17,ADVANCE ='no')    DREAL( v_x( i_x, i_y, i_z ) ), DIMAG( v_x( i_x, i_y, i_z ) )
      WRITE(73,f_c32p17,ADVANCE ='no')    DREAL( v_y( i_x, i_y, i_z ) ), DIMAG( v_y( i_x, i_y, i_z ) )
      WRITE(73,f_c32p17,ADVANCE ='yes')   DREAL( v_z( i_x, i_y, i_z ) ), DIMAG( v_z( i_x, i_y, i_z ) )

    END DO
    END DO
    END DO

    CLOSE(73)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real velocity. To read and plot velocity functions
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_3D )) // 'velocity_' &
              //TRIM( ADJUSTL( N_char ) ) // '_t_' // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where energy vs time will be written. With additional data

    OPEN( unit = 74, file = file_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO i_x = 0 , N_x - 1
    DO i_y = 0 , N_y - 1
    DO i_z = 0 , N_z - 1

      WRITE(74,f_d32p17,ADVANCE ='NO')    u_x( i_x, i_y, i_z )
      WRITE(74,f_d32p17,ADVANCE ='NO')    u_y( i_x, i_y, i_z )
      WRITE(74,f_d32p17,ADVANCE ='YES')   u_z( i_x, i_y, i_z )

    END DO
    END DO
    END DO

    CLOSE(74)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_running_status
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print the running status of the program, when called during debug
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN

      WRITE(*,'(A60)') TRIM( ADJUSTL( '-----------------------------------------------------------' ) )
      WRITE(*,'(A60)') TRIM( ADJUSTL( 'TIME  |   ENERGY   |   ENSTROPHY    |   INCOMPRESSIBILITY '  ) )
      WRITE(*,'(A60)') TRIM( ADJUSTL( '-----------------------------------------------------------' ) )

    END IF

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(F6.3,A3,F8.4,A3,F12.4,A7,E12.4)') time_now,'   ',energy,'   ',enstrophy,'     ',k_dot_v_norm
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF ( t_step .EQ. t_step_total ) THEN

      PRINT*,'-----------------------------------------------------------'

    END IF

  END

  SUBROUTINE print_error_timestep()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error if time step chosen is large
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40)')       TRIM( ADJUSTL( 'ERROR: TIME STEP TOO LARGE') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A40,F10.6)') TRIM( ADJUSTL( ' RESET THE TIME STEP (AT MAX) AS :') ),dt_max
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_error_nan()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error when NaN is encountered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40,F8.4)')  TRIM( ADJUSTL( 'ERROR: NAN ENCOUNTERED BEFORE T = ') ), time_now
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_error_incomp()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error when NaN is encountered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40,F8.4)')  TRIM( ADJUSTL( 'ERROR: INCOMPRESSIBILITY LOST BEFORE T = ') ), time_now
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_basicoutput
