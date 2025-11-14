module diff_op_mod
    use real_type_mod, only: dp
    use read_parameters_mod
    implicit none


    contains 


    subroutine linear_reconstruct(Nvar_local, u,u_plus,u_minus)
        implicit none
        integer, intent(in) :: Nvar_local
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar_local)
        real(dp), intent(out) :: u_plus(1-Ngz:Nx+Ngz,Nvar_local), u_minus(1-Ngz:Nx+Ngz,Nvar_local)
        real(dp) :: ratio(1-Ngz:Nx+Ngz,Nvar_local), phi(1-Ngz:Nx+Ngz,Nvar_local)


        ! Initialise arrays to zero
        ratio = 0.0_dp
        phi = 0.0_dp
        u_plus = 0.0_dp
        u_minus = 0.0_dp

        ! Compute ratio of slopes for limiter
        ! NB: We add small number to the denominator to 
        ! stop NaN issues when neighbouring points are close to equal
        ratio(0:Nx+1,:) = (u(0:Nx+1,:) - u(-1:Nx,:))/(u(1:Nx+2,:) - u(0:Nx+1,:) + 1.0e-16_dp) 

        ! Minmod limiter
        phi = max(0.0_dp,min(1.0_dp,ratio))

        ! MC Limiter
        ! phi = max(0.0_dp,min(min(2.0_dp*ratio,0.5_dp*(1.0_dp+ratio)),2.0_dp))

        ! Compute the slope-limited reconstruction
        u_minus(-1:Nx+1,:) = u(-1:Nx+1,:) + 0.5_dp*phi(-1:Nx+1,:)*(u(0:Nx+2,:) - u(-1:Nx+1,:))

        u_plus(-1:Nx,:) = u(0:Nx+1,:) - 0.5_dp*phi(0:Nx+1,:)*(u(1:Nx+2,:) - u(0:Nx+1,:))

    end subroutine linear_reconstruct

    subroutine local_lax_friedrichs(Nvar_local,flux_plus,flux_minus,cons_plus,cons_minus,char_speed,flux)
        implicit none
        integer, intent(in) :: Nvar_local
        real(dp), intent(in) :: flux_plus(1-Ngz:Nx+Ngz,Nvar_local), flux_minus(1-Ngz:Nx+Ngz,Nvar_local)
        real(dp), intent(in) :: cons_plus(1-Ngz:Nx+Ngz,Nvar_local), cons_minus(1-Ngz:Nx+Ngz,Nvar_local)
        real(dp), intent(in) :: char_speed(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: flux(1-Ngz:Nx+Ngz,Nvar_local)
        integer :: i
        

        !### NB: The sign of the dissipation changes depending on whether we evolving towards the future or the past!
        !### To the future = minus sign
        !### To the past = plus sign
        do i = 1, Nvar_local
            ! flux(:,i) = 0.5_dp*( (flux_plus(:,i) + flux_minus(:,i))  - char_speed(:,1)*(cons_plus(:,i)  - cons_minus(:,i))) ! Future flux

            flux(:,i) = 0.5_dp*( (flux_plus(:,i) + flux_minus(:,i))  + char_speed(:,1)*(cons_plus(:,i)  - cons_minus(:,i))) ! Past flux
        end do

    end subroutine local_lax_friedrichs


end module diff_op_mod