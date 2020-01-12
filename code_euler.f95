module parameters
! All the parameters for the simulation space are declared here
	implicit none
	double precision::lth,vol,vcy
	integer (kind=4)::ln2_N,N,no_t_step,N_half
	double complex::iota
	double precision(parameter)::two_pi
	
	contains 
	subroutine init_space_parameters
		implicit none
		ln2_N=4
		N=2**(ln2_N)
		N3=N**3
		! No of collocation points in physical space.
		N_half=N/2
		! No of collocation points in Fourier space
	end subroutine init_space_parameters
    subroutine init_time_parameters
        implicit none
        dt=0.0001
        no_t_step=10**6
        time_total=DFLOAT(no_t_step)*dt
    end subroutine init_time_parameters
end module parameters
