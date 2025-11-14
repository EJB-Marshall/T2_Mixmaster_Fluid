! Module containing the finite difference operations
module finite_difference_mod
    use real_type_mod, only: dp
    use read_parameters_mod
    implicit none

    contains

    ! subroutine FD_C4(dx,Nx,Ngz,u,du)
    !     ! Inputs:            
    !     ! dx - The spatial step size, real(dp)
    !     ! Nx - The number of grid points, integer
    !     ! Ngz - The number of ghost points, integer
    !     ! u - The function to be differentiated, real(dp) dimension(1,1-Ngz:Nx+Ngz)
    !     ! Output:
    !     ! du - The derivative, real(dp) dimension(1,1-Ngz:Nx+Ngz)
    !     implicit none
    !     real(dp), intent(in) :: dx 
    !     integer, intent(in) :: Nx, Ngz 
    !     real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,1)
    !     real(dp), intent(out) :: du(1-Ngz:Nx+Ngz,1)
    !     ! integer :: i

    !     ! do i = 1, Nx
    !     !     du(1,i) = (u(1,i-2) - 8.0_dp*u(1,i-1) + 8.0_dp*u(1,i+1) - u(1, i+2))/(12.0_dp*dx)
    !     ! end do

    !     du(1:Nx,:) = (u(-1:Nx-2,:) - 8.0_dp*u(0:Nx-1,:) + 8.0_dp*u(2:Nx+1,:) - u(3:Nx+2,:))/(12.0_dp*dx)

    
    ! end subroutine FD_C4



    function FD_C4(u) result(output)
        implicit none
        real(dp) :: u(1-Ngz:Nx+Ngz,1)
        real(dp) :: output(1-Ngz:Nx+Ngz,1)

        output(:,:) = 0.0_dp
        output(1:Nx,1) = (u(-1:Nx-2,1) - 8.0_dp*u(0:Nx-1,1) + 8.0_dp*u(2:Nx+1,1) - u(3:Nx+2,1))/(12.0_dp*dx)

    
    end function FD_C4

end module finite_difference_mod