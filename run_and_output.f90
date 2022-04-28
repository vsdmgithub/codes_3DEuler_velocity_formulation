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
! MODULE: run_and_output
! LAST MODIFIED: 29 JUNE 2020
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! RUN AND OUTPUT MODULE TO RUN, STUDY, PRINT OUTPUT FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE run_and_output
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module runs the Euler code, using solver module and other modules.
! 1. Main subroutine is 'time_forward_march' which moves through time steps
! 2. Subroutine 'spectral_data' collects spectral data and writes a shell wise energy.
! 3. Subroutine 'vortex_stretching' analyzes the way vorticity is stretched at every grid.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    USE solver
    USE paraview
    USE statistics
    ! List of modules refered and used.
    
    IMPLICIT NONE
    TYPE(VTR_file_handle)::fd_en          ! creating dataypes to store paraview files for 3d viewing   
    CONTAINS
    SUBROUTINE time_forward_march
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Loop of time steps, where at each step the spectral velocities
    ! are updated through any of the algoritm. Meanwhile, analysis and
    ! outputs are printed respectively.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        
        ALLOCATE(u_x(0:N-1,0:N-1,0:N-1),u_y(0:N-1,0:N-1,0:N-1),u_z(0:N-1,0:N-1,0:N-1))
        ALLOCATE(v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL init_initcondn
        ! Calls the subroutine to get a initial condition with norm_factor=1
        
        CALL compute_spectral_energy
        ! Gets the energy from spectral space
        
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  UNCOMMENT UNCOMMENT UNCOMMENT UNCOMMENT UNCOMMENT UNCOMMENT
        !-- N -- O -- R -- M -- A -- L -- I -- Z -- A -- T -- I -- O -- N --   
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
         
!        norm_factor=DSQRT(en_initial/norm_k)
!        ! Normalizing the norm_factor, so that we get energy='en_initial'
        
!        CALL init_initcondn
!        ! Calls it again to get normalized velocities.

!         CALL compute_spectral_energy
!        ! Gets the energy from spectral space   
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

        CALL write_parameters
        ! Writes the parameters

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        OPEN(unit=4,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name_en))//'.dat')
        ! File where energy vs time will be written. With additional data
        
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        !     E E U U L L E E R R    E E V V O O L L U U T T I I O O N N
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !             S        T         A         R       T
        ! 8888888888888888888888888888888888888888888888888888888888888888
 
        PRINT*,'-----------------------------------------------------------'
        PRINT*,'TIME  |  ENERGY  |  INCOMPRESSIBILITY  |  HELICITY'
        PRINT*,'-----------------------------------------------------------'
        DO t_step=0,t_step_total

            CALL step_to_time_convert(t_step,time_now)
            ! Converts the 't_step' to actual time 'time_now'

            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  A  N  A  L  Y  S  I  S       C   A   L   C  .
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            CALL spectral_data  ! ALSO PRINTS AT SAVE TIMES
            
            IF (viscosity_status .EQ. 0) THEN  
                CALL find_k_th
            END IF

            IF ((vortex_calc_status .EQ. 1) .AND. (MOD(t_step,t_step_save) .EQ. 0)) THEN  
                CALL strain_and_vorticity ! ALSO PRINTS AT SAVE TIMES
            END IF
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  E  N  E  R  G  Y      V  S      T  I  M  E
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                WRITE(4,f_d8p4,advance='no'),time_now
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp)
                WRITE(4,f_d32p17,advance='no'),Z_total
            IF (viscosity_status .EQ. 0) THEN  
                WRITE(4,f_d32p17,advance='no'),en_th
                WRITE(4,f_i4,advance='no'),k_th
            ELSE
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(k_kol:))
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(:k_kol-1))
            END IF

            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  U  R  G  I  N  G         T  I  M  E
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ((purging_status .EQ. 1) .AND. (MOD(t_step,t_step_purg) .EQ. 0)) then
            ! Its time to purge !
                v_x=purger*v_x
                v_y=purger*v_y
                v_z=purger*v_z
                purging_count=purging_count+1
            END IF
            IF (purging_status .EQ. 1 ) THEN
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(k_P:))
                WRITE(4,f_d32p17,advance='no'),flux_k_P(1)
                WRITE(4,f_d32p17,advance='no'),flux_k_P(2)
            END IF
            WRITE(4,f_i4,advance='yes'),purging_count
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ALLOCATE(v_x_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(conv_v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),conv_v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),conv_v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(conv_u_x(0:N-1,0:N-1,0:N-1),conv_u_y(0:N-1,0:N-1,0:N-1),conv_u_z(0:N-1,0:N-1,0:N-1))
            ALLOCATE(grad_x(0:N-1,0:N-1,0:N-1),grad_y(0:N-1,0:N-1,0:N-1),grad_z(0:N-1,0:N-1,0:N-1))
            ALLOCATE(dv1_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv1_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv1_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(v_x_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            CALL evolve_algorithm_rk4
            ! Updates v_x,v_y,v_z for next time step
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  D  E  A  L  L  O  C  A  T  I  O  N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            DEALLOCATE(v_x_dot,v_y_dot,v_z_dot)
            DEALLOCATE(conv_v_x,conv_v_y,conv_v_z)
            DEALLOCATE(conv_u_x,conv_u_y,conv_u_z)
            DEALLOCATE(grad_x,grad_y,grad_z)
            DEALLOCATE(dv1_x,dv2_x)
            DEALLOCATE(dv3_x,dv4_x)
            DEALLOCATE(dv1_y,dv2_y)
            DEALLOCATE(dv3_y,dv4_y)
            DEALLOCATE(dv1_z,dv2_z)
            DEALLOCATE(dv3_z,dv4_z)
            DEALLOCATE(v_x_temp,v_y_temp,v_z_temp)
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  D  E  B  U  G             F  O  R          N  a   N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF (MOD(t_step,t_step_debug) .EQ. 0) then
                CALL compute_spectral_energy
                CALL compute_real_energy
                CALL compute_helicity
                WRITE(*,'(F6.3,A3,F8.4,A3,F12.8,A8,F12.8)'),time_now,'   ',norm_k,'   ',k_dot_v,'         ',hely
                CALL compute_compressibility
                IF (count_nan .NE. 0) then
                    PRINT*,"NaN encountered before t=",time_now
                    exit
                    ! IF any NaN is encountered, the loop is exited without any further continuation.
               END IF
            END IF

        END DO
        PRINT*,'-----------------------------------------------------------'
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        !     E E U U L L E E R R    E E V V O O L L U U T T I I O O N N
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !                    E     N     D 
        ! 8888888888888888888888888888888888888888888888888888888888888888
        IF (pvd_status .EQ. 1) then
            CALL VTR_collect_file(FD=fd_en)
        END IF
        ! Standard routine to collect the pvd file
        
        CLOSE(4)
        
        state_sim=1
        ! Stating that the simulation has ended.
	END
    SUBROUTINE spectral_data
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This calculates energy, spectral shell wise. It goes through each
    ! spectral mode and puts the energy in the corresponding shell.
    ! This gives the ENERGY SPECTRUM.  When the time is right, it saves
    ! them too.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        Z_total=0.0D0
        en_sp=0.0D0
        ! Reset the array
        
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  E  N  E  R  G  Y     S  P  E  C  T  R  U  M     C  A   L   C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO j_z=-k_G,k_G
        DO j_y=-k_G,k_G
        DO j_x=1,k_G
            CALL mode_energy
            ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        END DO
        END DO
        j_x=0
        DO j_z=-k_G,-1
        DO j_y=-k_G,k_G
            CALL mode_energy             ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        END DO
        j_z=0
        DO j_y=0,k_G
            CALL mode_energy             ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  F  L  U  X       A  T      P  U  R  G  I  N  G      B  D  Y
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (purging_status .EQ. 1) THEN
             IF (t_step .EQ. 0) THEN
                  flux_k_P(1)=SUM(en_sp(:k_P))
             END IF
             flux_k_P(2)=(flux_k_P(1)-SUM(en_sp(:k_P)))/dt
             flux_k_P(1)=SUM(en_sp(:k_P))
        END IF
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  S  H  E  L  L      A  V  E  R  A  G  I  N  G  
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        en_sp_av(1)=hf*hf*(en_sp(2)+thr*en_sp(1))
        en_sp_av(max_shell_no)=hf*hf*(en_sp(max_shell_no-1)+thr*en_sp(max_shell_no))
        DO s=2,max_shell_no-1
            en_sp_av(s)=hf*hf*(en_sp(s-1)+en_sp(s+1)+two*en_sp(s))
        END DO
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  E  N  E  R  G  Y        O  U  T  P  U  T
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (MOD(t_step,t_step_save) .EQ. 0) then
            CALL write_spectral_energy
            ! Writes the energy spectrum.
        END IF
    END
    SUBROUTINE mode_energy
    ! CALL this to get the energy of that particular mode
        IMPLICIT NONE
        en_mode=CDABS(v_x(j_x,j_y,j_z))**two+CDABS(v_y(j_x,j_y,j_z))**two+CDABS(v_z(j_x,j_y,j_z))**two
    END
    SUBROUTINE find_k_th
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! NOTE: Only for Truncated/Purged systems.
    ! As thermalization begins, the bottom of the trench in the
    ! energy spectrum is denoted as the thermalization wavenumber.
    ! Trend is that starts near k_G and decreases all the way to
    ! the inertial range, which is a completely thermalized state.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        IF (purging_status .EQ. 1) THEN
            k_max=k_P-1
        ELSE
            k_max=k_G-1
        END IF
        
        IF (en_sp(k_max) .GT. tol) THEN
            act_therm=1
        END IF ! Activation of last shell

        IF (act_therm .EQ. 1) THEN
            k_th=MINLOC(en_sp_av(k_int:k_max-10), DIM=1)
            en_th=SUM(en_sp_av(k_th:k_max))
            ! Allowing the 'k_th' to reach till the integral scale 'k_int'
        ELSE
            k_th=k_max
            en_th=0.0D0
        END IF
    END
    SUBROUTINE strain_and_vorticity
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This calculates the vorticity in real space, and then strain tensor
    ! in real space. Then computes the dw/dt vector, and its component
    ! parallel to vorticity.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  A  L  L  O  C  A  T  I  O  N     F  O  R    C  A  L  C  . . 
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ALLOCATE(w_ux(0:N-1,0:N-1,0:N-1),w_uy(0:N-1,0:N-1,0:N-1),w_uz(0:N-1,0:N-1,0:N-1))
        ALLOCATE(str_xx(0:N-1,0:N-1,0:N-1),str_yy(0:N-1,0:N-1,0:N-1),str_zz(0:N-1,0:N-1,0:N-1))
        ALLOCATE(str_zx(0:N-1,0:N-1,0:N-1),str_xy(0:N-1,0:N-1,0:N-1),str_yz(0:N-1,0:N-1,0:N-1))
        ALLOCATE(w_mod(0:N-1,0:N-1,0:N-1),w_dot(0:N-1,0:N-1,0:N-1),w_dot_parallel(0:N-1,0:N-1,0:N-1))
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  V  O  R  T  I  C  I  T  Y       C  A  L  C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL fft_c2r(i*(k_y*v_z-k_z*v_y),i*(k_z*v_x-k_x*v_z),i*(k_x*v_y-k_y*v_x),N,Nh,w_ux,w_uy,w_uz)
        ! Spectral to Real Vorticity
                
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL fft_c2r(i*k_x*v_x,hf*i*(k_y*v_x+k_x*v_y),i*k_z*v_z,N,Nh,str_xx,str_xy,str_zz)
        CALL fft_c2r(i*k_y*v_y,hf*i*(k_y*v_z+k_z*v_y),hf*i*(k_x*v_z+k_z*v_x),N,Nh,str_yy,str_yz,str_zx)
               
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !   Q  U  A  T  E  R  N  I  O  N          C  A  L  C . .  .    .
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        w_mod=DSQRT(w_ux**two+w_uy**two+w_uz**two)
        ! Matrix consisisting of | \vec{w}| at grids.
        w_dot=DSQRT((str_xx*w_ux+str_xy*w_uy+str_zx*w_uz)**two+(str_xy*w_ux+str_yy*w_uy+str_yz*w_uz)**two &
        +(str_zx*w_ux+str_yz*w_uy+str_zz*w_uz)**two)
        ! Magnitude of 'S.w'
        w_dot_parallel=(str_xx*w_ux*w_ux+str_yy*w_uy*w_uy+str_zz*w_uz*w_uz+two*(str_xy*w_ux*w_uy+str_yz* &
        w_uy*w_uz+str_zx*w_ux*w_uz))/w_mod
        ! Projection of 'S.w'(pr_w_dot) along 'w'
        
!        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        !  D  E  A  L  L  O  C  A  T  I  O  N 
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DEALLOCATE(w_ux,w_uy,w_uz)
        DEALLOCATE(str_xx,str_yy,str_zz)
        DEALLOCATE(str_zx,str_xy,str_yz)
        ALLOCATE(quaternion(0:N-1,0:N-1,0:N-1))
        ALLOCATE(pdf_data(0:pdf_points),pdf_data_w(0:pdf_points),pdf_bins(0:pdf_points))
       
!        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!          P  .  D  .   F            O   U   T   P   U   T
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        quaternion_tag=1

        ! COSINE OF ALIGNMENT ANGLE
        quaternion=w_dot_parallel/w_dot

        min_val=-1.0D0
        max_val=1.0D0
        CALL pdf_function_adv(N,quaternion,w_mod,pdf_points,min_val,max_val,pdf_data,pdf_data_w,pdf_bins)
        CALL write_w_pdf
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        quaternion_tag=2

        ! VORTEX STRETCHING RATIO
        quaternion=w_dot_parallel/w_mod
        min_val=-20.0D0
        max_val=40.0D0
        CALL pdf_function_adv(N,quaternion,w_mod,pdf_points,min_val,max_val,pdf_data,pdf_data_w,pdf_bins)
        CALL write_w_pdf
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        quaternion_tag=3

        ! TWISTING RATIO
        quaternion=DSQRT(w_dot**two-w_dot_parallel**two)/w_mod
        min_val=0.0D0
        max_val=40.0D0
        CALL pdf_function_adv(N,quaternion,w_mod,pdf_points,min_val,max_val,pdf_data,pdf_data_w,pdf_bins)
        CALL write_w_pdf

!        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        !  D  E  A  L  L  O  C  A  T  I  O  N
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DEALLOCATE(w_mod,w_dot,w_dot_parallel)
        DEALLOCATE(quaternion)
        DEALLOCATE(pdf_data,pdf_data_w,pdf_bins)
    END
    SUBROUTINE write_w_pdf
    ! CALL this to WRITE the pdf of voriticity magnitude w.r.t the angle that w_dot and w makes
        IMPLICIT NONE
        CALL which_pdf
        file_name='pdf_'
        WRITE (file_time,f_d8p4),time_now
        ! Writes 'time_now' as a CHARACTER
        OPEN(unit=13,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name))//&
        TRIM(ADJUSTL(pdf_chr))//TRIM(ADJUSTL(file_time))//'.dat')
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  P  R  I  N   T          O  U  T  P  U  T   -  P.D.F
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO bin=1,pdf_points
             WRITE(13,f_d16p8,advance='no'),pdf_bins(bin)
             WRITE(13,f_d32p17,advance='no'),pdf_data(bin)
             WRITE(13,f_d32p17,advance='yes'),pdf_data_w(bin)
        END DO
        CLOSE(13)
    END
    SUBROUTINE which_pdf
    ! This determines which strain is handled, and which pdf is being calculated.
    IMPLICIT NONE
        IF (quaternion_tag .EQ. 1) THEN
            pdf_chr='proj_t_'
        ELSEIF (quaternion_tag .EQ. 2) THEN
            pdf_chr='alpha_t_'
        ELSE
            pdf_chr='chi_t_'
        END IF
    END
    SUBROUTINE write_spectral_energy
    ! CALL this to WRITE the spectral energy into a file with time index
        IMPLICIT NONE
        file_name='spectral_energy_t_'
        WRITE (file_time,f_d8p4),time_now
        ! Writes 'time_now' as a CHARACTER
        OPEN(unit=1,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name))//TRIM(ADJUSTL(file_time))//'.dat')
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  P  R  I  N   T          O  U  T  P  U  T   -   SPECTRUM
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO s1=1,k_G-1
             WRITE(1,f_d12p2,advance='no'),DFLOAT(s1)
             WRITE(1,f_d32p17,advance='no'),en_sp(s1)
             WRITE(1,f_d32p17,advance='yes'),en_sp_av(s1)
        END DO
        CLOSE(1)
    END
    END MODULE run_and_output
