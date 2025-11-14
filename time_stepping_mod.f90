! Module containing the time stepping routines
module time_stepping_mod
    use real_type_mod, only: dp
    use T2_fluid_mod, only: compute_rhs, CS_source, CS_source2
    use boundary_conditions_mod, only: periodic_bc
    use read_parameters_mod, only: Nx, Ngz, Nvar



    contains


    !#######################################
    !# `Classic' RK4 Solver
    !#######################################
    
    subroutine rk4(t0,y0,dt)
        implicit none
        real(dp), intent(in) :: t0,dt ! Current time, timestep 
        real(dp), intent(inout) :: y0(1-Ngz:Nx+Ngz,Nvar) ! Current values
        real(dp) :: k1(1-Ngz:Nx+Ngz,Nvar), k2(1-Ngz:Nx+Ngz,Nvar), k3(1-Ngz:Nx+Ngz,Nvar), k4(1-Ngz:Nx+Ngz,Nvar) ! Runge-Kutta Substeps
        real(dp) :: y1(1-Ngz:Nx+Ngz,Nvar), y2(1-Ngz:Nx+Ngz,Nvar), y3(1-Ngz:Nx+Ngz,Nvar) 
        real(dp) :: nu_norm, nu1, nu2, nu3, factor
        real(dp) :: eps
        integer :: i

        ! Step 1
        call compute_rhs(t0,y0,k1)
        y1 = y0 + 0.5*dt*k1


        ! Check for violations in fluid velocity constraint
        ! If necessary correct by re-scaling fluid components to max prescribed norm = 1.0 - eps
        eps = 1.0e-13_dp
        do i = 1, Nx
            nu1 = y1(i,10)
            nu2 = y1(i,11)
            nu3 = y1(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y1(i,10) = nu1 * factor
                y1(i,11) = nu2 * factor
                y1(i,12) = nu3 * factor
            end if
        end do
        
        call periodic_bc(Nx,Ngz,Nvar,y1) ! Update Boundary Conditions
        
        ! Step 2
        call compute_rhs(t0+0.5_dp*dt,y1,k2)
        y2 = y0 + 0.5*dt*k2

        do i = 1, Nx
            nu1 = y2(i,10)
            nu2 = y2(i,11)
            nu3 = y2(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y2(i,10) = nu1 * factor
                y2(i,11) = nu2 * factor
                y2(i,12) = nu3 * factor
            end if
        end do

        call periodic_bc(Nx,Ngz,Nvar,y2)

        ! Step 3
        call compute_rhs(t0+0.5_dp*dt,y2,k3)
        y3 = y0 + dt*k3

        do i = 1, Nx
            nu1 = y3(i,10)
            nu2 = y3(i,11)
            nu3 = y3(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y3(i,10) = nu1 * factor
                y3(i,11) = nu2 * factor
                y3(i,12) = nu3 * factor
            end if
        end do

        call periodic_bc(Nx,Ngz,Nvar,y3)

        ! Step 4
        call compute_rhs(t0+dt,y3,k4)

        
        ! Update y
        y0 = y0 + dt*((1.0_dp/6.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (1.0_dp/3.0_dp)*k3  + (1.0_dp/6.0_dp)*k4)

        do i = 1, Nx
            nu1 = y0(i,10)
            nu2 = y0(i,11)
            nu3 = y0(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y0(i,10) = nu1 * factor
                y0(i,11) = nu2 * factor
                y0(i,12) = nu3 * factor
            end if
        end do

        call periodic_bc(Nx,Ngz,Nvar,y0)


    end subroutine rk4




    !#######################################
    !# SSP RK3 Solver
    !#######################################

    subroutine rk3(t0,y0,dt)
        implicit none
        real(dp), intent(in) :: t0 ! Current time
        real(dp), intent(inout) :: dt ! Timestep
        real(dp), intent(inout) :: y0(1-Ngz:Nx+Ngz,Nvar) ! Current values
        real(dp) :: k1(1-Ngz:Nx+Ngz,Nvar), k2(1-Ngz:Nx+Ngz,Nvar), k3(1-Ngz:Nx+Ngz,Nvar) ! Runge-Kutta Substeps
        real(dp) :: y1(1-Ngz:Nx+Ngz,Nvar), y2(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: nu_norm, nu1, nu2, nu3, factor
        real(dp) :: eps
        integer :: i

        ! Step 1
        call compute_rhs(t0,y0,k1)
        y1 = y0 + dt*k1
        ! dt = max(dt,-CS_source*0.05_dp) ! Norm preserving timestep needs current value of time derivative!


        eps = 1.0e-13_dp
        do i = 1, Nx
            nu1 = y1(i,10)
            nu2 = y1(i,11)
            nu3 = y1(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y1(i,10) = nu1 * factor
                y1(i,11) = nu2 * factor
                y1(i,12) = nu3 * factor
            end if
        end do

        call periodic_bc(Nx,Ngz,Nvar,y1) ! Update Boundary Conditions

         
        
        ! Step 2
        call compute_rhs(t0 + dt,y1,k2)
        y2 = y0 + 0.25_dp*dt*(k1 + k2)


        do i = 1, Nx
            nu1 = y2(i,10)
            nu2 = y2(i,11)
            nu3 = y2(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y2(i,10) = nu1 * factor
                y2(i,11) = nu2 * factor
                y2(i,12) = nu3 * factor
            end if
        end do

        call periodic_bc(Nx,Ngz,Nvar,y2)

        

        ! Step 3
        call compute_rhs(t0 + 0.5_dp*dt,y2,k3)


        ! Update y
        y0 = y0 + dt*((1.0_dp/6.0_dp)*k1 + (1.0_dp/6.0_dp)*k2 + (2.0_dp/3.0_dp)*k3)


        do i = 1, Nx
            nu1 = y0(i,10)
            nu2 = y0(i,11)
            nu3 = y0(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y0(i,10) = nu1 * factor
                y0(i,11) = nu2 * factor
                y0(i,12) = nu3 * factor
            end if
        end do


        call periodic_bc(Nx,Ngz,Nvar,y0)



    end subroutine rk3


    !#######################################
    !# RK1 (Forward Euler)
    !#######################################

    subroutine rk1(t0,y0,dt)
        implicit none
        real(dp), intent(in) :: t0 ! Current time 
        real(dp), intent(inout) :: dt ! Timestep
        real(dp), intent(inout) :: y0(1-Ngz:Nx+Ngz,Nvar) ! Current values
        real(dp) :: k1(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: y1(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: nu_norm, nu1, nu2, nu3, factor
        real(dp) :: eps
        integer :: i

        ! Step 1
        call compute_rhs(t0,y0,k1)
        ! print *, "CS_source is", CS_source
        ! dt = max(dt,-CS_source*0.015_dp) ! Norm preserving timestep needs current value of time derivative!

        y1 = y0 + dt*k1

        eps = 1.0e-13_dp
        do i = 1, Nx
            nu1 = y1(i,10)
            nu2 = y1(i,11)
            nu3 = y1(i,12)

            nu_norm = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)

            if (nu_norm > 1.0_dp - eps) then

                factor = (1.0_dp - eps) / nu_norm

                y1(i,10) = nu1 * factor
                y1(i,11) = nu2 * factor
                y1(i,12) = nu3 * factor
            end if
        end do

        print *, maxval(sqrt(y1(1:Nx,10)**(2.0_dp) + y1(1:Nx,11)**(2.0_dp)  + y1(1:Nx,12)**(2.0_dp)))

        y0 = y1

        call periodic_bc(Nx,Ngz,Nvar,y0) ! Update Boundary Conditions



    end subroutine rk1



end module time_stepping_mod