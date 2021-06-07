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
! MODULE: system_advectionsolver
! LAST MODIFIED: 29 JUNE 2020
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D EULER EQUATIONS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advectionsolver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and moves it one step forward in euler equation using the given algorithm, uses FFTW.
! This has RK4, RK2, AB2, AB4 algorithm. Additionally subroutines 'compute_spectral_energy', 'compute_real_energy', 'compute_compressibility'
! gives energy in spectral, real space and gives the output of \vec{k}\cdot \vec{v}.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_variables
	USE system_fftw

	IMPLICIT NONE
	! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::adv_u_x,adv_u_y,adv_u_z
  ! Real advection term matrix

  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::grad_x,grad_y,grad_z
  ! Calculates gradient in real space from FFT

	! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::adv_v_x,adv_v_y,adv_v_z
	! Spectral advection term matrix

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

	SUBROUTINE allocate_advectionsolver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE(adv_u_x(0:N-1,0:N-1,0:N-1),adv_u_y(0:N-1,0:N-1,0:N-1),adv_u_z(0:N-1,0:N-1,0:N-1))
		ALLOCATE(grad_x(0:N-1,0:N-1,0:N-1),grad_y(0:N-1,0:N-1,0:N-1),grad_z(0:N-1,0:N-1,0:N-1))
		ALLOCATE(adv_v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),adv_v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),adv_v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

	END

	SUBROUTINE allocate_advectionsolver_RK4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for RK4 algoritm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(dv1_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv1_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv1_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(v_x_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

	END

	SUBROUTINE advectionsolver_RK4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
		v_x_temp = v_x
		v_y_temp = v_y
		v_z_temp = v_z
		CALL time_increment_RK1() ! This call provides \vec{dv} for the existing \vec{v}
		v_x      = v_x_temp + hf * dv1_x
		v_y      = v_y_temp + hf * dv1_y
		v_z      = v_z_temp + hf * dv1_z
		CALL time_increment_RK2()
		v_x      = v_x_temp + hf * dv2_x
		v_y      = v_y_temp + hf * dv2_y
		v_z      = v_z_temp + hf * dv2_z
		CALL time_increment_RK3()
		v_x      = v_x_temp + dv3_x
		v_y      = v_y_temp + dv3_y
		v_z      = v_z_temp + dv3_z
		CALL time_increment_RK4()
		! Final increment for 'v(k)'
		v_x      = v_x_temp + ( dv1_x + two * dv2_x + two * dv3_x + dv4_x ) / six
		v_y      = v_y_temp + ( dv1_y + two * dv2_y + two * dv3_y + dv4_y ) / six
		v_z      = v_z_temp + ( dv1_z + two * dv2_z + two * dv3_z + dv4_z ) / six

	END

	SUBROUTINE allocate_advectionsolver_AB4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for AB4 algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE(v_x_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_x_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_x_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_x_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_y_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_z_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_x_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

	END

	SUBROUTINE advectionsolver_AB4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE AB4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Adam Bashforth 4th Order algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		IF ( t_step .EQ. 0 ) THEN

			! Using initial condition to get -3 step value
			CALL time_derivative_AB()
			v_x_dot_m3 = v_x_dot
			v_y_dot_m3 = v_y_dot
			v_z_dot_m3 = v_z_dot

			CALL allocate_advectionsolver_RK4
			! Allocating RK4 for first 3 steps

			CALL advectionsolver_RK4_algorithm

		ELSE IF ( t_step .EQ. 1 ) THEN

			! Filling up for -2 step value
			CALL time_derivative_AB()
			v_x_dot_m2 = v_x_dot
			v_y_dot_m2 = v_y_dot
			v_z_dot_m2 = v_z_dot

			CALL advectionsolver_RK4_algorithm

		ELSE IF ( t_step .EQ. 2 ) THEN

			! Filling up for -1 step value
			CALL time_derivative_AB()
			v_x_dot_m1 = v_x_dot
			v_y_dot_m1 = v_y_dot
			v_z_dot_m1 = v_z_dot

			CALL advectionsolver_RK4_algorithm

			CALL deallocate_advectionsolver_RK4
			! No need for RK4 anymore

		ELSE

			CALL time_derivative_AB()

			v_x_pred = v_x + dt * ( - 9.0D0 * v_x_dot_m3 + 37.0D0 * v_x_dot_m2 - 59.0D0 * v_x_dot_m1 + 55.0D0 * v_x_dot ) / 24.0D0
			v_y_pred = v_y + dt * ( - 9.0D0 * v_y_dot_m3 + 37.0D0 * v_y_dot_m2 - 59.0D0 * v_y_dot_m1 + 55.0D0 * v_y_dot ) / 24.0D0
			v_z_pred = v_z + dt * ( - 9.0D0 * v_z_dot_m3 + 37.0D0 * v_z_dot_m2 - 59.0D0 * v_z_dot_m1 + 55.0D0 * v_z_dot ) / 24.0D0

			CALL time_derivative_AB_pred()

			v_x      = v_x + dt * ( v_x_dot_m2 - 5.0D0 * v_x_dot_m1 + 19.0D0 * v_x_dot + 9.0D0 * v_x_dot_m3 ) / 24.0D0
			v_y      = v_y + dt * ( v_y_dot_m2 - 5.0D0 * v_y_dot_m1 + 19.0D0 * v_y_dot + 9.0D0 * v_y_dot_m3 ) / 24.0D0
			v_z      = v_z + dt * ( v_z_dot_m2 - 5.0D0 * v_z_dot_m1 + 19.0D0 * v_z_dot + 9.0D0 * v_z_dot_m3 ) / 24.0D0
			! Predicted 'v_dot' is stored in 'v_dot_m3' - to save space :)

			! Shiting the known velocities for next step
			v_x_dot_m3 = v_x_dot_m2
			v_x_dot_m2 = v_x_dot_m1
			v_x_dot_m1 = v_x_dot

			v_y_dot_m3 = v_y_dot_m2
			v_y_dot_m2 = v_y_dot_m1
			v_y_dot_m1 = v_y_dot

			v_z_dot_m3 = v_z_dot_m2
			v_z_dot_m2 = v_z_dot_m1
			v_z_dot_m1 = v_z_dot

		END IF

	END

	SUBROUTINE time_derivative_AB()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First FFT spectral to real velocity
		CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )

		! u_x gradient in real space
		CALL fft_c2r( i * k_x * v_x, i * k_y * v_x, i * k_z * v_x, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in x direction
		adv_u_x = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_y, i * k_y * v_y, i * k_z * v_y, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_y = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_z, i * k_y * v_z, i * k_z * v_z, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_z = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c( adv_u_x, adv_u_y, adv_u_z, N, Nh, adv_v_x, adv_v_y, adv_v_z )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		v_x_dot = - truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		v_y_dot = - truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		v_z_dot = - truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )

	END

	SUBROUTINE time_derivative_AB_pred()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First FFT spectral to real velocity
		CALL fft_c2r( v_x_pred, v_y_pred, v_z_pred, N, Nh, u_x, u_y, u_z )

		! u_x gradient in real space
		CALL fft_c2r( i * k_x * v_x_pred, i * k_y * v_x_pred, i * k_z * v_x_pred, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in x direction
		adv_u_x = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_y_pred, i * k_y * v_y_pred, i * k_z * v_y_pred, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_y = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_z_pred, i * k_y * v_z_pred, i * k_z * v_z_pred, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_z = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c( adv_u_x, adv_u_y, adv_u_z, N, Nh, adv_v_x, adv_v_y, adv_v_z )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		v_x_dot_m3 = - truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		v_y_dot_m3 = - truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		v_z_dot_m3 = - truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )

	END

	SUBROUTINE time_increment_RK1()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL advection
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.
		dv1_x = - dt * truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		dv1_y = - dt * truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		dv1_z = - dt * truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )
	END

	SUBROUTINE time_increment_RK2()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL advection
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.
		dv2_x = - dt * truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		dv2_y = - dt * truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		dv2_z = - dt * truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )
	END

	SUBROUTINE time_increment_RK3()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL advection
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.
		dv3_x = - dt * truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		dv3_y = - dt * truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		dv3_z = - dt * truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )
	END

	SUBROUTINE time_increment_RK4()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL advection
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.
		dv4_x = - dt * truncator * ( proj_xx * adv_v_x + proj_xy * adv_v_y + proj_zx * adv_v_z )
		dv4_y = - dt * truncator * ( proj_xy * adv_v_x + proj_yy * adv_v_y + proj_yz * adv_v_z )
		dv4_z = - dt * truncator * ( proj_zx * adv_v_x + proj_yz * adv_v_y + proj_zz * adv_v_z )
	END

	SUBROUTINE advection
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to give spectral advection term using v_k.
	! 1. First v 			--> u  FFT is done
	! 2. Next i*k*v 	--> du/dx  FFT is done
	! 3. Next u.du/dx --> Fourier(u.du/dx)
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First FFT spectral to real velocity
		CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )

		! u_x gradient in real space
		CALL fft_c2r( i * k_x * v_x, i * k_y * v_x, i * k_z * v_x, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in x direction
		adv_u_x = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_y, i * k_y * v_y, i * k_z * v_y, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_y = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! u_y gradient in real space
		CALL fft_c2r( i * k_x * v_z, i * k_y * v_z, i * k_z * v_z, N, Nh, grad_x, grad_y, grad_z )

		! u.Nabla(u) term in z direction
		adv_u_z = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c( adv_u_x, adv_u_y, adv_u_z, N, Nh, adv_v_x, adv_v_y, adv_v_z )

  END

	SUBROUTINE deallocate_advectionsolver_RK4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(dv1_x,dv2_x)
		DEALLOCATE(dv3_x,dv4_x)
		DEALLOCATE(dv1_y,dv2_y)
		DEALLOCATE(dv3_y,dv4_y)
		DEALLOCATE(dv1_z,dv2_z)
		DEALLOCATE(dv3_z,dv4_z)
		DEALLOCATE(v_x_temp,v_y_temp,v_z_temp)

	END

	SUBROUTINE deallocate_advectionsolver_AB4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(v_x_pred,v_y_pred,v_z_pred)
		DEALLOCATE(v_x_dot_m1,v_x_dot_m2,v_x_dot_m3)
		DEALLOCATE(v_y_dot_m1,v_y_dot_m2,v_y_dot_m3)
		DEALLOCATE(v_z_dot_m1,v_z_dot_m2,v_z_dot_m3)

	END

	SUBROUTINE deallocate_advectionsolver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(adv_v_x,adv_v_y,adv_v_z)
		DEALLOCATE(adv_u_x,adv_u_y,adv_u_z)
		DEALLOCATE(grad_x,grad_y,grad_z)

	END

 END MODULE system_advectionsolver
