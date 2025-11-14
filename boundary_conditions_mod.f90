! Module containing the boundary conditions
module boundary_conditions_mod
    use real_type_mod, only: dp
    implicit none

    contains 

    subroutine periodic_bc(Nx,Ngz,Nvar,u)
        ! Inputs:            
        ! dx - The spatial step size, real(dp)
        ! Nx - The number of grid points, integer
        ! Nvar - The number of variables we are updating
        ! u - The function we apply the BCs to, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        ! Output:
        ! u - The function with updated boundaries, real(dp) dimension(1,1-Ngz:Nx+Ngz)

        implicit none
        integer, intent(in) :: Nx, Ngz, Nvar
        real(dp), intent(inout) :: u(1-Ngz:Nx+Ngz,Nvar)


        ! Right Boundary
        u(Nx+1:Nx+Ngz,:) = u(1:Ngz,:)

        ! Left Boundary
        u(1-Ngz:0,:) = u(Nx-Ngz+1:Nx,:)



    end subroutine periodic_bc

end module boundary_conditions_mod