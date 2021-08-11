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
! MODULE: system_advvariables
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED VARIABLES TO DO ANALYSIS IN 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_basicfunctions

  IMPLICIT NONE
  ! _________________________
  ! VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  :: k_P
  INTEGER(KIND=4)  :: gr1_size,gr2_size
  INTEGER(KIND=4)  :: ind_1,ind_2
  DOUBLE PRECISION :: dum_double
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::shell_en_matrix, shell_es_matrix
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::shell_en_XS1, shell_es_XS1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::shell_en_XS2, shell_es_XS2
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE  ::sh_gr1_x,sh_gr1_y,sh_gr1_z
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE  ::sh_gr2_x,sh_gr2_y,sh_gr2_z
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::k_x_axis,k_y_axis,k_z_axis

  CONTAINS

  SUBROUTINE allocate_shell_grid
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::Delta_k, k_mod

    Delta_k = one
    ind_1   = 0
    ind_2   = 0
    k_P     = k_G - 5

    DO j_x = kMin_x, kMax_x
  	DO j_y = kMin_y, kMax_y
  	DO j_z = kMin_z, kMax_z

      k_mod    = DSQRT( k_2( j_x, j_y, j_z ) )

      IF ( DABS( k_mod - k_G ) .LT. Delta_k ) THEN
        ind_1  = ind_1 + 1
      END IF

      IF ( DABS( k_mod - k_P ) .LT. Delta_k ) THEN
        ind_2  = ind_2 + 1
      END IF

    END DO
    END DO
    END DO

    gr1_size = ind_1
    gr2_size = ind_2

    ind_1   = 0
    ind_2   = 0
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( sh_gr1_x( gr1_size ),sh_gr1_y( gr1_size ),sh_gr1_z( gr1_size )   )
    ALLOCATE( sh_gr2_x( gr2_size ),sh_gr2_y( gr2_size ),sh_gr2_z( gr2_size )   )

    DO j_x = kMin_x, kMax_x
  	DO j_y = kMin_y, kMax_y
  	DO j_z = kMin_z, kMax_z

      k_mod                    = DSQRT( k_2( j_x, j_y, j_z ) )

      IF ( DABS( k_mod - k_G ) .LT. Delta_k ) THEN
        ind_1             = ind_1 + 1
        sh_gr1_x( ind_1 ) = j_x
        sh_gr1_y( ind_1 ) = j_y
        sh_gr1_z( ind_1 ) = j_z
      END IF

      IF ( DABS( k_mod - k_P ) .LT. Delta_k ) THEN
        ind_2             = ind_2 + 1
        sh_gr2_x( ind_2 ) = j_x
        sh_gr2_y( ind_2 ) = j_y
        sh_gr2_z( ind_2 ) = j_z
      END IF

    END DO
    END DO
    END DO

    ALLOCATE( shell_en_XS1( gr1_size ) )
    ALLOCATE( shell_es_XS1( gr1_size ) )
    ALLOCATE( shell_en_XS2( gr2_size ) )
    ALLOCATE( shell_es_XS2( gr2_size ) )

    ALLOCATE( k_x_axis( kMin_x:kMax_x ) )
    ALLOCATE( k_y_axis( kMin_y:kMax_y ) )
    ALLOCATE( k_z_axis( kMin_z:kMax_z ) )

    DO j_x = kMin_x, kMax_x
      k_x_axis( j_x ) = K_scale_x * DBLE( j_x )
    END DO

    DO j_y = kMin_y, kMax_y
      k_y_axis( j_y ) = K_scale_y * DBLE( j_y )
    END DO

    DO j_z = kMin_z, kMax_z
      k_z_axis( j_z ) = K_scale_z * DBLE( j_z )
    END DO

    ALLOCATE( shell_en_matrix( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( shell_es_matrix( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

  END

  SUBROUTINE deallocate_shell_grid
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to deallocate
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D  E  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( sh_gr1_x, sh_gr1_y, sh_gr1_z )
    DEALLOCATE( sh_gr2_x, sh_gr2_y, sh_gr2_z )
    DEALLOCATE( shell_en_XS1 )
    DEALLOCATE( shell_en_XS2 )
    DEALLOCATE( shell_es_XS1 )
    DEALLOCATE( shell_es_XS2 )
    DEALLOCATE( k_x_axis, k_y_axis, k_z_axis )
    DEALLOCATE( shell_en_matrix )
    DEALLOCATE( shell_es_matrix )

  END

END MODULE system_advvariables
