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
  DOUBLE PRECISION :: loc_stretching, vx_stretching_max, loc_vx_stretching_max
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::str_xx,str_yy,str_zz,str_xy,str_yz,str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_xx,bck_str_yy,bck_str_zz,bck_str_xy,bck_str_yz,bck_str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_opr,w_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_vx_stretching,vx_stretching

  CONTAINS

  SUBROUTINE allocate_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(w_mod_2(0:N-1,0:N-1,0:N-1))
    ALLOCATE(str_xx(0:N-1,0:N-1,0:N-1),str_yy(0:N-1,0:N-1,0:N-1),str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(str_zx(0:N-1,0:N-1,0:N-1),str_xy(0:N-1,0:N-1,0:N-1),str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(vx_stretching(0:N-1,0:N-1,0:N-1))

  END

  SUBROUTINE allocate_bck_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the bckground strain operator array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::arg, r_local

    r_local = 3.5D0 * dx
    ! Size of the region to integrate and consider as local

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(bck_str_opr(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(bck_str_xx(0:N-1,0:N-1,0:N-1),bck_str_yy(0:N-1,0:N-1,0:N-1),bck_str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_str_zx(0:N-1,0:N-1,0:N-1),bck_str_xy(0:N-1,0:N-1,0:N-1),bck_str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_vx_stretching(0:N-1,0:N-1,0:N-1))

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      arg                          = r_local * DSQRT( k_2( i_x, i_y, i_z ) )

      bck_str_opr( i_x, i_y, i_z ) = thr * ( DSIN( arg ) - arg * DCOS( arg ) ) /  (arg ** thr)

    END DO
    END DO
    END DO

    bck_str_opr( 0, 0, 0 ) = one
    ! Zero mode which is background

  END

	SUBROUTINE deallocate_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to strain
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(w_mod_2)
		DEALLOCATE(str_xx,str_yy,str_zz)
		DEALLOCATE(str_xy,str_yz,str_zx)
		DEALLOCATE(vx_stretching)

	END

	SUBROUTINE deallocate_bck_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to background strain
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(bck_str_xx,bck_str_yy,bck_str_zz)
    DEALLOCATE(bck_str_xy,bck_str_yz,bck_str_zx)
		DEALLOCATE(bck_vx_stretching)
		DEALLOCATE(bck_str_opr)

	END

END MODULE system_advvariables
