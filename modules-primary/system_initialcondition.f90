  ! <f
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
! MODULE: system_initialcondition
! LAST MODIFIED: 21 June 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR 3D EULER EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
  ! </f>

MODULE system_initialcondition
  ! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Different initial conditions are coded here
! 1. Exponentially decaying spectrum
! 2. Thermalized spectrum
! 3. Kolmogorov like spectrum
! 4. Taylor-Green Initial Condition
! 5. Kida Peltz Initial Condition
! 6. ABC Flow initial condition
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! </f>
  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_fftw_adv
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT  NONE

  CONTAINS

  SUBROUTINE init_initcondn
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for velocity, Choose one.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    CALL init_fft_size( N_x, N_y, N_z )
    ! Initializing the size of domain for FFT's in simulation - One time procedure.
    ! REF-> <<< system_fftw_adv >>>

    ! Initializing the initial velocity (spectral) and projecting it so that the flow is incompressible.

    ! CALL IC_exp_decaying_spectrum(energy_initial)
    ! Generic randomized initial condition, with energy mainly in integral scale (spectrally)

    ! CALL IC_Kolmogorov_spectrum(energy_initial)
    ! Generic initial condition, with energy mainly in inertial range with a k^-(5/3) spectrum.

    ! CALL IC_perfect_thermalized_spectrum(energy_initial)
    ! Create its own thermalized spectrum by equiparition, (no permanence of large eddies in this case)

    ! CALL IC_vortex_sheet(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from IC_exp_decaying_spectrum

    ! CALL IC_vortex_sheet_compression(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a local disturbance

    CALL IC_disturbedsheet(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a  local disturbance

    ! CALL IC_vortex_sheet_with_TG(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from Taylor Green

    ! CALL IC_vortex_tube(energy_initial)
    ! Creates a vortex tube at z = 0, along z direction.
    ! With a background field from IC_exp_decaying_spectrum

    ! CALL IC_vortex_cylinder(energy_initial)
    ! Creates a vortex cylinder, from Bell and Marcus paper 1990
    ! With a background field from Taylor Green

    ! CALL IC_from_file_spectral
    ! Read from file.
    ! *****Check whether file is available already.

    ! CALL IC_from_file_real
    ! Read from file.
    ! *****Check whether file is available already.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! UNCOMMENT TO TRUNCATE IF NEEDED - (most of I.C are already truncated)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    v_x = truncator * v_x
    v_y = truncator * v_y
    v_z = truncator * v_z

  END
  ! </f>

  SUBROUTINE IC_exp_decaying_spectrum(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition IN which only few large modes are excited.
  ! CALL this to give 3 components of COMPLEX spectral velocity at spectral grid 'i_x,i_y,i_z'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::integral_exponent

    IC_type = 'EXP-DEC'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    k_integral        = 2
    ! Integral scale wavenumber

    integral_exponent = 2
    ! The power in the spectrum E(k) for k<k_integral
    ! Generally either 2 or 4 or 6. But has to be a even number

    CALL normalization_exponential_spectrum( integral_exponent, k_integral, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_integral )

      V_k_mod = norm_factor * norm_const * k_ratio**( hf * integral_exponent - 1 ) &
                * DEXP( - qtr * integral_exponent * ( k_ratio ** two ) )

      V_k(1)  = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)  = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)  = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is Incompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given Initial velocity
      ! But those planes require special attention, becoz, their Inversion lies IN the same plane,
      ! so conjugates must be placed accordingly.

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
        v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
      END IF

      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)

        v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
        v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END
  ! </f>

  SUBROUTINE IC_Kolmogorov_spectrum(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition with kolmogorov spectrum model, referred from Pope's Turbulence.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE PRECISION             ::factor_integral,factor_dissipation
    DOUBLE PRECISION             ::eleven_by_twelve
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::integral_exponent
    INTEGER(KIND=4)              ::k_kol

    IC_type = 'KOL-INE'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    k_integral        = FLOOR( DLOG( DBLE( N_max ) ) / DLOG( 4.0D0 ) ) - 1
    ! Integral scale wavenumber

    k_kol             = FLOOR ( DBLE( N_max ) / 4.0D0 ) - 1
    ! End of kolmogorov spectrum

    integral_exponent = 2
    ! The power in the spectrum E(k) for k<k_integral
    ! Generally either 2 or 4 or 6. But has to be a even number

    eleven_by_twelve  = 11.0D0 /  12.0D0

    CALL normalization_kolmogorov_spectrum( k_integral, k_kol, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi                = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta              = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph                 = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_integral )
      factor_integral    = one
      factor_dissipation = one

      IF ( k_ratio .LT. 1 ) THEN
        CALL kolmogorov_spectrum_integralscale_subpart(k_ratio,integral_exponent,factor_integral)
        ! REF-> <<< system_auxilaries >>>
      END IF

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_kol )

      IF ( k_ratio .GT. qtr ) THEN
        CALL kolmogorov_spectrum_dissipationscale_subpart(k_ratio,factor_dissipation)
        ! REF-> <<< system_auxilaries >>>
      END IF

      V_k_mod            = norm_factor * norm_const * ( k_2( i_x, i_y, i_z ) ** ( - eleven_by_twelve ) ) &
                          * factor_integral * factor_dissipation

      V_k(1)             = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)             = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)             = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given INitial velocity
      ! But those planes require special attention, becoz, their INversion lies IN the same plane,
      ! so conjugates must be placed accordINgly.

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
        v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
      END IF

      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)

        v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
        v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END
  ! </f>

  SUBROUTINE IC_perfect_thermalized_spectrum(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for thermalized spectrum
  ! where equipartition is given.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::phi,theta,norm_const
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE COMPLEX,DIMENSION(3)::V_k

    IC_type = 'THER-K2'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    norm_const = DSQRT( thr / ( two_pi * DBLE( k_G ** 3 ) ) )

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      V_k(1)  = norm_factor * norm_const * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)  = norm_factor * norm_const * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)  = norm_factor * norm_const * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given INitial velocity
      ! But those planes require special attention, becoz, their INversion lies IN the same plane,
      ! so conjugates must be placed accordINgly.

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
       v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                  proj_zx( i_x, i_y, i_z ) * V_k(3)
       v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                  proj_yz( i_x, i_y, i_z ) * V_k(3)
       v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                  proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
       v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
       v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
       v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
      END IF

      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
       v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_zx( i_x, i_y, i_z) * V_k(3)
       v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_yz( i_x, i_y, i_z) * V_k(3)
       v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_zz( i_x, i_y, i_z) * V_k(3)

       v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
       v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
       v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
       v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_zx( i_x, i_y, i_z) * V_k(3)
       v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_yz( i_x, i_y, i_z) * V_k(3)
       v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                  proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END
  ! </f>

  SUBROUTINE IC_vortex_sheet(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a background field.
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    DOUBLE PRECISION::energy_sheet,energy_ratio
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.2D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio = 0.001D0
    ! Percentage of energy in Background field

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = u0 * ( two + DTANH( - c_factor * DBLE( i_x - ( N_x / 4 ) ) ) &
      + DTANH( c_factor * DBLE( i_x - 3 * ( N_x / 4 ) ) ) )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting the real velocity to add the sheet

    u_y = u_y + u_sheet_y

    CALL fft_r2c( u_y, v_y )
    ! FFT spectral to real velocity

    IC_type = 'VOR-SHT'

    DEALLOCATE(u_sheet_y)

  END
  ! </f>

  SUBROUTINE IC_vortex_sheet_compression(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a simple strain compressing the sheet at the center
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    DOUBLE PRECISION::energy_sheet,energy_ratio,energy_comp
    DOUBLE PRECISION::gaus1,gaus3,x_gr1,x_gr3,y_gr,z_gr,k_beta
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.25D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio = 0.005D0
    ! Percentage of energy in compression field

    i_x0 = INT( N_x / 2 )
    i_y0 = INT( N_y / 2 )
    i_z0 = INT( N_z / 2)
    i_x1 = 1 * INT( N_x / 4 )
    i_x3 = 3 * INT( N_x / 4 )

    k_beta = 4.0D0
    ! Spread of the gaussian damping

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = u0 * ( two + DTANH( - c_factor * DBLE( i_x - 1 * ( N_x / 4 ) ) ) &
                                              + DTANH( + c_factor * DBLE( i_x - 3 * ( N_x / 4 ) ) ) )

      x_gr1                      = DBLE( i_x - i_x1 ) * dx
      x_gr3                      = DBLE( i_x - i_x3 ) * dx
      y_gr                       = DBLE( i_y - i_y0 ) * dy
      z_gr                       = DBLE( i_z - i_z0 ) * dz

      ! gaus1                      = DEXP( - hf * ( k_beta ** two ) * ( x_gr1 ** two + y_gr ** two + z_gr ** two ) )
      ! gaus3                      = DEXP( - hf * ( k_beta ** two ) * ( x_gr3 ** two + y_gr ** two + z_gr ** two ) )

      gaus1                      = DEXP( - hf * ( k_beta ** two ) * ( x_gr1 ** two + z_gr ** two ) )
      gaus3                      = DEXP( - hf * ( k_beta ** two ) * ( x_gr3 ** two + z_gr ** two ) )

      u_x( i_x, i_y, i_z )       = - x_gr1 * gaus1 + x_gr3 * gaus3
      u_y( i_x, i_y, i_z )       = zero
      u_z( i_x, i_y, i_z )       = + z_gr  * gaus1 - z_gr  * gaus3

    END DO
    END DO
    END DO

    CALL fft_r2c( u_sheet_y, v_y )
    ! FFT spectral to real velocity

    CALL fourier_smoothing
    ! Smooths the spectrum near the truncation to remove any gibbsian oscillations

    CALL fft_c2r( v_y, u_sheet_y )
    ! Real velocity to spectral velocity

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    CALL compute_projected_velocity
    ! Makes the compressing field divergenceless

    v_x( 0, 0, 0 ) = c0
    v_y( 0, 0, 0 ) = c0
    v_z( 0, 0, 0 ) = c0

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting real velocity

    energy_comp  = hf * SUM( ( u_x ** two ) + ( u_y ** two ) + ( u_z ** two ) ) / N3
    u0           = DSQRT( energy_ratio * energy_input / energy_comp )
    u_x          = u0 * u_x
    u_y          = u0 * u_y
    u_z          = u0 * u_z
    ! Normalisation of comp flow

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    u_y          = u_y + u_sheet_y
    ! ! Combining both flows (superposition)

    CALL fft_r2c( u_y, v_y )
    ! FFT spectral to real velocity

    IC_type      = 'SHT_CMP'

  END
  ! </f>

  SUBROUTINE IC_disturbedsheet(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a small localised disturbance, parametrized by
  ! an angle \Psi, along which we expect oscillations to arise
  ! Amount of energy in the disturbance is adjustable by 'energy_ratio'
  ! Spread of the disturbance is characterised by the wavenumber 'k_beta' - essentially in the large scale
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    DOUBLE PRECISION::psi,cs,sn
    DOUBLE PRECISION::x_ro1,x_ro3,y_ro1,y_ro3,z_ro1,z_ro3,y_ro,z_ro
    DOUBLE PRECISION::ux_ro1,ux_ro3,uy_ro1,uy_ro3,uz_ro1,uz_ro3
    DOUBLE PRECISION::energy_sheet,energy_ratio,energy_dist
    DOUBLE PRECISION::gaus1,gaus3,x_gr1,x_gr3,y_gr,z_gr,k_beta
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0                 = one
    ! Normalizing parameter

    smooth_pm          = 0.200D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor           = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio       = 0.005D0
    ! Percentage of energy in compression field

    psi                = 30.0D0 * ( two_pi / 360.0D0 )
    cs                 = DCOS( psi )
    sn                 = DSIN( psi )
    ! Angle at which the compression is oriented

    i_x0               = INT( N_x / 2 )
    i_y0               = INT( N_y / 2 )
    i_z0               = INT( N_z / 2)
    i_x1               = 1 * INT( N_x / 4 )
    i_x3               = 3 * INT( N_x / 4 )

    k_beta             = 4.0D0
    ! Spread of the gaussian damping

    ! Disturbance in XY Plane
    ! DO i_z = 0, N_z - 1
  	! DO i_y = 0, N_y - 1
  	! DO i_x = 0, N_x - 1
    !
    !   u_sheet_y( i_x, i_y, i_z ) = u0 * ( two + DTANH( - c_factor * DBLE( i_x - 1 * ( N_x / 4 ) ) ) &
    !                                           + DTANH( + c_factor * DBLE( i_x - 3 * ( N_x / 4 ) ) ) )
    !
    !   x_gr1                      = DBLE( i_x - i_x1 ) * dx
    !   x_gr3                      = DBLE( i_x - i_x3 ) * dx
    !   y_gr                       = DBLE( i_y - i_y0 ) * dy
    !   z_gr                       = DBLE( i_z - i_z0 ) * dz
    !
    !   z_ro                       = z_gr
    !
    !   x_ro1                      = x_gr1 * cs + y_gr * sn
    !   y_ro1                      = y_gr  * cs - x_gr1 * sn
    !   gaus1                      = DEXP( - hf * ( k_beta ** two ) * ( x_ro1 ** two + y_ro1 ** two + z_ro ** two ) )
    !   ux_ro1                     = - x_ro1 * gaus1
    !   uy_ro1                     = + y_ro1 * gaus1
    !
    !   x_ro3                      = x_gr3 * cs + y_gr * sn
    !   y_ro3                      = y_gr  * cs - x_gr3 * sn
    !   gaus3                      = DEXP( - hf * ( k_beta ** two ) * ( x_ro3 ** two + y_ro3 ** two + z_ro ** two ) )
    !   ux_ro3                     = + x_ro3 * gaus3
    !   uy_ro3                     = - y_ro3 * gaus3
    !
    !   u_x( i_x, i_y, i_z )       = ( ux_ro1 + ux_ro3 ) * cs - ( uy_ro1 + uy_ro3 ) * sn
    !   u_y( i_x, i_y, i_z )       = ( ux_ro1 + ux_ro3 ) * sn + ( uy_ro1 + uy_ro3 ) * cs
    !   u_z( i_x, i_y, i_z )       = zero
    !
    ! END DO
    ! END DO
    ! END DO
    !
    ! Disturbance in XZ Plane
    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = u0 * ( one + DTANH( - c_factor * DBLE( i_x - 1 * ( N_x / 4 ) ) ) &
                                              + DTANH( + c_factor * DBLE( i_x - 3 * ( N_x / 4 ) ) ) )

      x_gr1                      = DBLE( i_x - i_x1 ) * dx
      x_gr3                      = DBLE( i_x - i_x3 ) * dx
      y_gr                       = DBLE( i_y - i_y0 ) * dy
      z_gr                       = DBLE( i_z - i_z0 ) * dz

      y_ro                       = y_gr

      x_ro1                      = x_gr1 * cs + z_gr * sn
      z_ro1                      = z_gr  * cs - x_gr1 * sn
      gaus1                      = DEXP( - hf * ( k_beta ** two ) * ( x_ro1 ** two + y_ro ** two + z_ro1 ** two ) )
      ux_ro1                     = - x_ro1 * gaus1
      uz_ro1                     = + z_ro1 * gaus1

      x_ro3                      = x_gr3 * cs + z_gr * sn
      z_ro3                      = z_gr  * cs - x_gr3 * sn
      gaus3                      = DEXP( - hf * ( k_beta ** two ) * ( x_ro3 ** two + y_ro ** two + z_ro3 ** two ) )
      ux_ro3                     = + x_ro3 * gaus3
      uz_ro3                     = - z_ro3 * gaus3

      u_x( i_x, i_y, i_z )       = ( ux_ro1 + ux_ro3 ) * cs - ( uz_ro1 + uz_ro3 ) * sn
      u_y( i_x, i_y, i_z )       = zero
      u_z( i_x, i_y, i_z )       = ( ux_ro1 + ux_ro3 ) * sn + ( uz_ro1 + uz_ro3 ) * cs

    END DO
    END DO
    END DO

    CALL fft_r2c( u_sheet_y, v_y )
    ! FFT spectral to real velocity

    CALL fourier_smoothing
    ! Smooths the spectrum near the truncation to remove any gibbsian oscillations

    CALL fft_c2r( v_y, u_sheet_y )
    ! Real velocity to spectral velocity

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    CALL compute_projected_velocity
    ! Makes the compressing field divergenceless

    v_x( 0, 0, 0 ) = c0
    v_y( 0, 0, 0 ) = c0
    v_z( 0, 0, 0 ) = c0

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting real velocity

    energy_dist  = hf * SUM( ( u_x ** two ) + ( u_y ** two ) + ( u_z ** two ) ) / N3
    u0           = DSQRT( energy_ratio * energy_input / energy_dist )
    u_x          = u0 * u_x
    u_y          = u0 * u_y
    u_z          = u0 * u_z
    ! Normalisation of comp flow

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    u_y          = u_y + u_sheet_y
    ! ! Combining both flows (superposition)

    CALL fft_r2c( u_y, v_y )
    ! FFT spectral to real velocity

    IC_type      = 'DIS_30'

  END
  ! </f>

  SUBROUTINE IC_vortex_sheet_with_TG(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a TAYLOR GREEN FLOW as background.
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3
    DOUBLE PRECISION::energy_sheet,energy_ratio,energy_TG
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.8D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio = 0.005D0
    ! Percentage of energy in Background field

    ! i_x0 =     INT( N_x /  8)
    i_x0 = 0
    i_y0 = 0
    i_z0 = 0
    i_x1 = 1 * INT( N_x / 4 )
    i_x3 = 3 * INT( N_x / 4 )

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = two + DTANH( - c_factor * DBLE( i_x - i_x1 ) ) &
                                       + DTANH( + c_factor * DBLE( i_x - i_x3 ) )

      u_x( i_x, i_y, i_z )       = + DCOS( DBLE( i_x - i_x0 ) * dx )&
                                   * DSIN( DBLE( i_y - i_y0 ) * dy )&
                                   * DCOS( DBLE( i_z - i_z0 ) * dz )
      u_y( i_x, i_y, i_z )       = - DSIN( DBLE( i_x - i_x0 ) * dx )&
                                   * DCOS( DBLE( i_y - i_y0 ) * dy )&
                                   * DCOS( DBLE( i_z - i_z0 ) * dz )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    energy_TG    = hf * SUM( ( u_x ** two ) + ( u_y ** two ) ) / N3
    u0           = DSQRT( energy_ratio * energy_input / energy_TG )
    u_x          = u0 * u_x
    u_y          = u0 * u_y
    ! Normalisation of TG flow

    u_y = u_y + u_sheet_y
    ! Combining both flows (superposition)

    IC_type = 'VOR-STG'
    ! vortex sheet with Taylor Green

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    DEALLOCATE(u_sheet_y)

  END
  ! </f>

  SUBROUTINE IC_vortex_tube(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a background field.
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm
    DOUBLE PRECISION::tube_y0,tube_z0
    DOUBLE PRECISION::tube_y,tube_z
    DOUBLE PRECISION::energy_tube,energy_ratio
    DOUBLE PRECISION::u_ang,radius,arg
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_tube_y,u_tube_z

    ALLOCATE(u_tube_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE(u_tube_z( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.04D0
    ! How thick the sheet is, smaller the parameter thicker it is

    energy_ratio = 0.01D0
    ! Percentage of energy in Background field

    tube_y0      = DBLE( N_y / 2 )
    tube_z0      = DBLE( N_z / 2 )
    ! Center of the tube

    DO i_y = 0, N_y - 1
    DO i_z = 0, N_z - 1

      tube_y                  = DBLE( i_y ) - tube_y0
      tube_z                  = DBLE( i_z ) - tube_z0
      radius                  = DSQRT( tube_y ** two + tube_z ** two )
      arg                     = radius * smooth_pm * two_pi / thr
      u_ang                   = arg * DEXP( - hf * ( arg ** two ) )
      u_tube_y( :, i_y, i_z ) = - u_ang * tube_z / radius
      u_tube_z( :, i_y, i_z ) = + u_ang * tube_y / radius

    END DO
    END DO

    energy_tube  = hf * SUM( u_tube_y ** two + u_tube_z ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_tube )
    u_tube_y     = u0 * u_tube_y
    u_tube_z     = u0 * u_tube_z
    ! Normalization of tube

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting the real velocity to add the tube velocity

    u_y = u_y + u_tube_y
    u_z = u_z + u_tube_z

    CALL fft_r2c( u_y, v_y )
    CALL fft_r2c( u_z, v_z )
    ! FFT spectral to real velocity

    CALL compute_projected_velocity
    ! Projects the velocity to remove some compressibility

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_initial / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

    IC_type = 'VOR-TUB'

    DEALLOCATE( u_tube_y, u_tube_z )

  END
  ! </f>

  SUBROUTINE IC_vortex_cylinder(energy_input)
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex cylinder
  ! Thickness of the cylinder and the strength of the perturbation is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    DOUBLE PRECISION::smooth_pm,c_factor
    DOUBLE PRECISION::r_cyl,r_ind,r_gri
    DOUBLE PRECISION::G_arg,pert,k_beta

    smooth_pm                            = 0.25D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor                             = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor                   => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    r_cyl                                = two_pi / 6.0D0
    ! Radius of the cylinder

    r_gri                                = r_cyl / dy
    ! Grid no for the radius of the cylinder

    pert                                 = 0.10
    ! Perturbation of the cylinder, pert =0 means no stationary flow

    k_beta                               = 2.0D0
    ! Width of the gaussian perturbation

    i_x0 = N_x / 2
    i_y0 = N_y / 2
    i_z0 = N_z / 2

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      r_ind                = DSQRT( DBLE( ( i_y - i_y0 ) ** two + ( i_z - i_z0 ) ** two ) )
      u_x( i_x, i_y, i_z ) = DTANH( c_factor * ( r_gri - r_ind ) )

      u_y( i_x, i_y, i_z ) = zero

      G_arg                = hf * ( k_beta ** two ) * ( ( ( i_x - i_x0 ) * dx ) ** two + ( ( i_y - i_y0 ) * dy ) ** two )
      u_z( i_x, i_y, i_z ) = pert * DEXP( - G_arg )

    END DO
    END DO
    END DO

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

    IC_type = 'VOR-CYL'

  END
  ! </f>

  SUBROUTINE IC_from_file_spectral
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in spectral space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::real_part,imag_part
    CHARACTER(LEN=80)::IC_file_name

    IC_type = 'INP-SPE'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = 'spectral_velocity_' // TRIM( ADJUSTL( N_char ) ) // '_input.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 43, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    DO i_x = kMin_x, kMax_x
    DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_x( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_y( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='YES') real_part, imag_part
      v_z( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )

    END DO
    END DO
    END DO

    CLOSE(43)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END
  ! </f>

  SUBROUTINE IC_from_file_real
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in real space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::IC_file_name

    IC_type = 'INP-REA'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = 'velocity_' // TRIM( ADJUSTL( N_char ) ) // '_input.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 44, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    DO i_x = 0 , N_x - 1
    DO i_y = 0 , N_y - 1
    DO i_z = 0 , N_z - 1

      READ(44,f_d32p17,ADVANCE = 'NO')  u_x( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'NO')  u_y( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'YES') u_z( i_x, i_y, i_z)

    END DO
    END DO
    END DO

    CLOSE(44)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! FFT spectral to real velocity

  END
  ! </f>

  SUBROUTINE compute_energy_spectral_data
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check the presence of NaN in your spectral velocity data (v(k)),
  ! and also the L2 norm or the Kinetic energy.
  ! NOTE: Count certain modes once, certain modes half (owing to 1/2 factor)
  ! in the first loop i_x=0 plane is left. later it is considered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    energy      = zero

    DO i_x      =  kMin_x + 1, kMax_x - 1
    DO i_y      =  kMin_y    , kMax_y
    DO i_z      =  kMin_z    , kMax_z
      energy    = energy + CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x         =   kMin_x
    DO i_y      =   kMin_y , kMax_y
    DO i_z      =   kMin_z , kMax_z
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

    i_x         =   kMax_x
    DO i_y      =   kMin_y , kMax_y
    DO i_z      =   kMin_z , kMax_z
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

  END
  ! </f>

  SUBROUTINE compute_projected_velocity
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to project the spectral velocity, to make it incompressible
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::v_P_x,v_P_y,v_P_z
    ALLOCATE( v_P_x( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_P_y( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_P_z( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

    v_P_x = v_x
    v_P_y = v_y
    v_P_z = v_z

    v_x   = proj_xx * v_P_x + proj_xy * v_P_y + proj_zx * v_P_z
    v_y   = proj_xy * v_P_x + proj_yy * v_P_y + proj_yz * v_P_z
    v_z   = proj_zx * v_P_x + proj_yz * v_P_y + proj_zz * v_P_z

    DEALLOCATE(v_P_x,v_P_y,v_P_z)

  END
  ! </f>

  SUBROUTINE fourier_smoothing
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to smoothen the spectral velocities to deprive of any Gibbsian oscillations
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::smoothing_exp,k_cutoff,arg

    smoothing_exp = 6.0D0
    ! Power in the exponential, atleast 2, can go upto 32

    k_cutoff = DBLE( N_min / 3 )
    ! cutoff wavenumber where the filter starts to act vigorously.

    DO j_x = kMin_x, kMax_x
  	DO j_y = kMin_y, kMax_y
  	DO j_z = kMin_z, kMax_z

      arg = DEXP( - ( DSQRT( k_2( j_x, j_y, j_z ) ) / k_cutoff ) ** smoothing_exp )

      v_x( j_x, j_y, j_z ) = v_x( j_x, j_y, j_z ) * arg
      v_y( j_x, j_y, j_z ) = v_y( j_x, j_y, j_z ) * arg
      v_z( j_x, j_y, j_z ) = v_z( j_x, j_y, j_z ) * arg

    END DO
    END DO
    END DO

  END
  ! </f>

  SUBROUTINE check_compressibility
  ! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility condition. Sums over all residues
  ! of incompressibility and prints it. Of order 10^(-12).
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_dot_v_norm   = zero

    DO i_x         = kMin_x + 1, kMax_x - 1
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                           k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                           k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x            = kMin_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                               k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                               k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    i_x            = kMax_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                               k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                               k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    k_dot_v_norm = DSQRT( k_dot_v_norm )

    print*,'COMPRESSIBILITY = ',k_dot_v_norm

  END
  ! </f>

END MODULE system_initialcondition
