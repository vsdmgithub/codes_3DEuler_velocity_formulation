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

! #########################
! MODULE: system_interactionsolver
! LAST MODIFIED: 29 JUNE 2020
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D EULER EQUATIONS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_interactionsolver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and moves it one step forward in euler equation using the given algorithm, uses FFTW.
! This has RK4, AB4 algorithm.
! The spectral equation is
! dv_i(k)/dt = P_ij . ( u(x) X w(x) )_j(k)
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicvariables
	USE system_fftw_adv

	IMPLICIT NONE
	! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::uXw_x,uXw_y,uXw_z
	! Real cross product between velocity and vorticity

	! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::vXw_x, vXw_y, vXw_z
  ! Spectral cross product between velocity and vorticity

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m1,v_y_dot_m1,v_z_dot_m1
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m2,v_y_dot_m2,v_z_dot_m2
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m3,v_y_dot_m3,v_z_dot_m3
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot,v_y_dot,v_z_dot
	! Spectral derivatives for AB4

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_pred,v_y_pred,v_z_pred
  ! Spectral velocity predictor for AB4 matrix

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv1_x,dv2_x,dv3_x,dv4_x,dv1_y,dv2_y
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv3_y,dv4_y,dv1_z,dv2_z,dv3_z,dv4_z
  ! Intermediate matrices for RK4 algorithm

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_temp,v_y_temp,v_z_temp
  ! temporary matrices to store velocities during RK4 algorithm


	CONTAINS

	SUBROUTINE allocate_interactionsolver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ALLOCATE( grad_1( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( grad_2( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

	END

	SUBROUTINE merger_TG_VX_interaction
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! INTERACTION TERM
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		w_ux_md = w_ux_md + dw_ux
		w_uy_md = w_uy_md + dw_uy
		w_uz_md = w_uz_md + dw_uz

	END

	SUBROUTINE interactionsolver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! INTERACTION TERM
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		CALL fft_c2r( i * k_y * w_ux, grad_1 )
		CALL fft_c2r( i * k_z * u_x, grad_2 )

		w_ux_md = w_ux_md + dt * (- u_y_VX * grad_1 + w_uz_VX * grad_2)

		CALL fft_c2r( i * k_y * w_uy, grad_1 )
		CALL fft_c2r( i * k_y * u_z, grad_2 )

		w_uy_md = w_uy_md + dt * (- u_y_VX * grad_1 + w_uz_VX * grad_2)

		CALL fft_c2r( i * k_y * w_uz, grad_1 )

		w_uz_md = w_uz_md + dt * (- u_y_VX * grad_1 - w_uz_VX_der * u_x)

	END

	SUBROUTINE deallocate_interactionsolver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		DEALLOCATE( grad_1, grad_2 )

	END

 END MODULE system_interactionsolver
