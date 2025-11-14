! Module containing compute_rhs subroutine for the T2 Einstein-Euler equations 
module T2_fluid_mod
    use real_type_mod, only: dp
    use read_parameters_mod
    use finite_difference_mod
    use diff_op_mod
    implicit none
    real(dp) :: CS_x, CS_source, CS_source2

    contains 


    !########################################################
    !# Compute Characteristic Speeds
    !########################################################

    subroutine compute_characteristic_speeds(prims_plus, prims_minus, E11_plus, E11_minus, char_speed)
        implicit none
        ! real(dp), intent(in) :: K
        real(dp), intent(in) :: prims_plus(1-Ngz:Nx+Ngz,4), prims_minus(1-Ngz:Nx+Ngz,4)
        real(dp), intent(in) :: E11_plus(1-Ngz:Nx+Ngz,1), E11_minus(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: char_speed(1-Ngz:Nx+Ngz,1)
        real(dp) :: nu1_plus(1-Ngz:Nx+Ngz,1), nu2_plus(1-Ngz:Nx+Ngz,1), nu3_plus(1-Ngz:Nx+Ngz,1)
        real(dp) :: nu1_minus(1-Ngz:Nx+Ngz,1), nu2_minus(1-Ngz:Nx+Ngz,1), nu3_minus(1-Ngz:Nx+Ngz,1)
        real(dp) :: lam1_plus(1-Ngz:Nx+Ngz,1), lam1_minus(1-Ngz:Nx+Ngz,1), lam2_plus(1-Ngz:Nx+Ngz,1), lam2_minus(1-Ngz:Nx+Ngz,1)
        real(dp) :: lam3_plus(1-Ngz:Nx+Ngz,1), lam3_minus(1-Ngz:Nx+Ngz,1), lam4_plus(1-Ngz:Nx+Ngz,1), lam4_minus(1-Ngz:Nx+Ngz,1)

        nu1_plus(:,1) = prims_plus(:,1)
        nu2_plus(:,1) = prims_plus(:,2)
        nu3_plus(:,1) = prims_plus(:,3)

        nu1_minus(:,1) = prims_minus(:,1)
        nu2_minus(:,1) = prims_minus(:,2)
        nu3_minus(:,1) = prims_minus(:,3)


        lam1_plus = E11_plus*nu1_plus !# \lambda_{1} has multiplicity 2
        lam1_minus = E11_minus*nu1_minus
        
        lam2_plus = E11_plus*((1.0_dp-K)*nu1_plus - sqrt(K*(1.0_dp-nu1_plus**(2.0_dp) - nu2_plus**(2.0_dp) - nu3_plus**(2.0_dp))*(1.0_dp-nu1_plus**(2.0_dp) &
                    - K*(nu2_plus**(2.0_dp) + nu3_plus**(2.0_dp)))))/(1.0_dp - K*(nu1_plus**(2.0_dp) + nu2_plus**(2.0_dp) + nu3_plus**(2.0_dp)))

        lam2_minus = E11_minus*((1.0_dp-K)*nu1_minus - sqrt(K*(1.0_dp-nu1_minus**(2.0_dp) - nu2_minus**(2.0_dp) - nu3_minus**(2.0_dp))*(1.0_dp-nu1_minus**(2.0_dp) &
                     - K*(nu2_minus**(2.0_dp) + nu3_minus**(2.0_dp)))))/(1.0_dp - K*(nu1_minus**(2.0_dp) + nu2_minus**(2.0_dp) + nu3_minus**(2.0_dp)))

        lam3_plus = E11_plus*((1.0_dp-K)*nu1_plus + sqrt(K*(1.0_dp-nu1_plus**(2.0_dp) - nu2_plus**(2.0_dp) - nu3_plus**(2.0_dp))*(1.0_dp-nu1_plus**(2.0_dp) &
                    - K*(nu2_plus**(2.0_dp) + nu3_plus**(2.0_dp)))))/(1.0_dp - K*(nu1_plus**(2.0_dp) + nu2_plus**(2.0_dp) + nu3_plus**(2.0_dp)))

        lam3_minus = E11_minus*((1.0_dp-K)*nu1_minus + sqrt(K*(1.0_dp-nu1_minus**(2.0_dp)-nu2_minus**(2.0_dp) - nu3_minus**(2.0_dp))*(1.0_dp-nu1_minus**(2.0_dp) &
                    - K*(nu2_minus**(2.0_dp) + nu3_minus**(2.0_dp)))))/(1.0_dp - K*(nu1_minus**(2.0_dp) + nu2_minus**(2.0_dp) + nu3_minus**(2.0_dp)))
        

        !# The eigenvalues for the gravitational equations are \pm E^{1}_{1}
        !# but we only need the absolute value, so just compute one of them
        lam4_plus = E11_plus 

        lam4_minus = E11_minus


        ! Compute the largest value of the characteristic speeds at each cell
        char_speed = max(abs(lam1_plus), abs(lam1_minus), &
                        abs(lam2_plus), abs(lam2_minus), & 
                        abs(lam3_plus), abs(lam3_minus), &
                        abs(lam4_plus), abs(lam4_minus))
        
        ! Compute the largest value of the characteristic speeds over the *entire* domain
        CS_x = maxval(char_speed)
        CS_x = max(CS_x, 1.0e-12_dp) ! Make sure CS_x doesn't return zero


    end subroutine compute_characteristic_speeds


    !########################################################################################################
    !# Compute Fluxes
    !########################################################################################################

    subroutine compute_flux_grav(cons_grav_plus,cons_grav_minus,flux_grav_plus, flux_grav_minus)
        implicit none
        real(dp), intent(in) :: cons_grav_plus(1-Ngz:Nx+Ngz,4), cons_grav_minus(1-Ngz:Nx+Ngz,4)
        real(dp), intent(out) :: flux_grav_plus(1-Ngz:Nx+Ngz,4), flux_grav_minus(1-Ngz:Nx+Ngz,4)
        real(dp) :: Sigma_Minus_plus(1-Ngz:Nx+Ngz,1), Sigma_Times_plus(1-Ngz:Nx+Ngz,1), N_Minus_plus(1-Ngz:Nx+Ngz,1), N_Times_plus(1-Ngz:Nx+Ngz,1)
        real(dp) :: Sigma_Minus_minus(1-Ngz:Nx+Ngz,1), Sigma_Times_minus(1-Ngz:Nx+Ngz,1), N_Minus_minus(1-Ngz:Nx+Ngz,1), N_Times_minus(1-Ngz:Nx+Ngz,1)

        !### Plus Flux
        Sigma_Minus_plus(:,1) = cons_grav_plus(:,1) 
        Sigma_Times_plus(:,1) = cons_grav_plus(:,2) 
        N_Minus_plus(:,1) = cons_grav_plus(:,3) 
        N_Times_plus(:,1) = cons_grav_plus(:,4) 

        flux_grav_plus(:,1) = N_Times_plus(:,1) !# Sigma_Minus flux
        flux_grav_plus(:,2) = -N_Minus_plus(:,1) !# Sigma_Times flux
        flux_grav_plus(:,3) = -Sigma_Times_plus(:,1) !# N_Minus flux
        flux_grav_plus(:,4) = Sigma_Minus_plus(:,1) !# N_Times flux


        !### Minus Flux
        Sigma_Minus_minus(:,1) = cons_grav_minus(:,1) 
        Sigma_Times_minus(:,1) = cons_grav_minus(:,2) 
        N_Minus_minus(:,1) = cons_grav_minus(:,3) 
        N_Times_minus(:,1) = cons_grav_minus(:,4) 


        flux_grav_minus(:,1) = N_Times_minus(:,1)
        flux_grav_minus(:,2) = -N_Minus_minus(:,1)
        flux_grav_minus(:,3) = -Sigma_Times_minus(:,1)
        flux_grav_minus(:,4) = Sigma_Minus_minus(:,1)


    end subroutine compute_flux_grav



    !########################################################################################################
    !# Compute Fluid Derivative Coefficients
    !########################################################################################################

    !### This function computes the coefficients of the non-conservative derivative terms in the fluid equations
    !### The equations are of the form \del_{t}U = B(U)\del_{x}U
    !### where U = (log(T^{00}), nu1, nu2, nu3)^{T}


    !### Thus, B(U) is a 4x4 matrix.
    !### B(1,1) is the coefficient of the \del_{x}log(T^{00}) term in the log(T^{00}) evo eqn
    !### B(1,2) is the coefficient of the \del_{x}log(nu1) term in the log(T^{00}) evo eqn,
    !### B(2,2) is the coefficient of the \del_{x}log(nu1) term in the nu1 evo eqn 
    !### etc. etc. 

    subroutine compute_deriv_coefficients(prims_fluid,E11,B_coeff)
        implicit none
        real(dp), intent(in) :: prims_fluid(1-Ngz:Nx+Ngz,4), E11(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: B_coeff(1-Ngz:Nx+Ngz,4,4)
        real(dp) :: nu1(1-Ngz:Nx+Ngz,1), nu2(1-Ngz:Nx+Ngz,1), nu3(1-Ngz:Nx+Ngz,1), nu_norm(1-Ngz:Nx+Ngz,1)

        nu1(:,1) = prims_fluid(:,1)
        nu2(:,1) = prims_fluid(:,2)
        nu3(:,1) = prims_fluid(:,3)

        nu_norm = nu1**(2.0_dp) + nu2**(2.0_dp)  + nu3**(2.0_dp)

        !# Initialise B matrix
        B_coeff(:,:,:) = 0.0_dp


        !### log(tau) Equation Coefficients

        B_coeff(:,1,1) = -(1.0_dp+K)/(1.0_dp+K*nu_norm(:,1))*nu1(:,1)*E11(:,1)

        B_coeff(:,1,2) = -((1.0_dp+K)*(1.0_dp+K*nu_norm(:,1)-2.0_dp*K*nu1(:,1)**(2.0_dp)))/((1.0_dp + K*nu_norm(:,1))**(2.0_dp))*E11(:,1)

        B_coeff(:,1,3) = (2.0_dp*K*(1.0_dp+K)*nu1(:,1)*nu2(:,1))/((1.0_dp + K*nu_norm(:,1))**(2.0_dp))*E11(:,1)

        B_coeff(:,1,4) = (2.0_dp*K*(1.0_dp + K)*nu1(:,1)*nu3(:,1))/((1.0_dp + K*nu_norm(:,1))**(2.0_dp))*E11(:,1)

        !### nu1 Equation Coefficients

        B_coeff(:,2,1) = (K*(1.0_dp + K*nu_norm(:,1)**(2.0_dp) + nu1(:,1)**(2.0_dp)*(-1.0_dp + K + 2.0_dp*(1.0_dp + K)*nu_norm(:,1)) &
                            -nu_norm(:,1)*(1.0_dp + K + (1.0_dp + 3.0_dp*K)*nu1(:,1)**(2.0_dp))))/((1.0_dp + K)*(-1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        B_coeff(:,2,2) = (nu1(:,1)*(1.0_dp - 3.0_dp*K +(1.0_dp + K)*K*nu_norm(:,1) + 2.0_dp*K*(1.0_dp - K)*nu1(:,1)**(2.0_dp)))&
                        /(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        B_coeff(:,2,3) = (2.0_dp*K*(-1.0_dp + K*nu_norm(:,1) - (K-1.0_dp)*nu1(:,1)**(2.0_dp))*nu2(:,1))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        B_coeff(:,2,4) = (2.0_dp*K*(-1.0_dp + K*nu_norm(:,1) - (K-1.0_dp)*nu1(:,1)**(2.0_dp))*nu3(:,1))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        !### nu2 Equation Coefficients

        B_coeff(:,3,1) = (K*nu1(:,1)*nu2(:,1)*(1.0_dp - K)*(nu_norm(:,1) -1.0_dp))/((1.0_dp+K)*(-1.0_dp+K*nu_norm(:,1)))*E11(:,1)

        B_coeff(:,3,2) = (K*nu2(:,1)*((nu_norm(:,1)-1.0_dp)*(1.0_dp + K*nu_norm(:,1)) - 2.0_dp*(K-1.0_dp)*nu1(:,1)**(2.0_dp)))&
                        /(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        B_coeff(:,3,3) = nu1(:,1)*(-1.0_dp - (2.0_dp*(K-1.0_dp)*K*nu2(:,1)**(2.0_dp))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp)))*E11(:,1)

        B_coeff(:,3,4) = (2.0_dp*(1.0_dp - K)*nu1(:,1)*nu2(:,1)*nu3(:,1))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        !### nu3 Equation Coefficients

        B_coeff(:,4,1) = (K*nu1(:,1)*nu3(:,1)*(1.0_dp - K)*(nu_norm(:,1) - 1.0_dp))/((1.0_dp + K)*(-1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        B_coeff(:,4,2) = (K*nu3(:,1)*((nu_norm(:,1)-1.0_dp)*(1.0_dp + K*nu_norm(:,1)) - 2.0_dp*(K-1.0_dp)*nu1(:,1)**(2.0_dp)))&
                        /(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        B_coeff(:,4,3) = (2.0_dp*(1.0_dp-K)*K*nu1(:,1)*nu2(:,1)*nu3(:,1))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp))*E11(:,1)

        B_coeff(:,4,4) = nu1(:,1)*( -1.0_dp - (2.0_dp*(K-1.0_dp)*K*nu3(:,1)**(2.0_dp))/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**(2.0_dp)))*E11(:,1)




        !#### Test by importing directly from Mathematica:

        !### log(tau) Equation Coefficients

        ! B_coeff(:,1,1) = -(((1.0_dp + K)*nu1(:,1))/(1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        ! B_coeff(:,1,2) = -(((1.0_dp + K)*(1.0_dp + K*nu_norm(:,1) - 2.0_dp*K*nu1(:,1)**2.0_dp))/(1.0_dp + K*nu_norm(:,1))**2.0_dp)*E11(:,1)

        ! B_coeff(:,1,3) = ((2.0_dp*K*(1.0_dp + K)*nu1(:,1)*nu2(:,1))/(1.0_dp + K*nu_norm(:,1))**2.0_dp)*E11(:,1)

        ! B_coeff(:,1,4) = ((2.0_dp*K*(1.0_dp + K)*nu1(:,1)*nu3(:,1))/(1.0_dp + K*nu_norm(:,1))**2.0_dp)*E11(:,1)

        ! !### nu1 Equation Coefficients

        ! B_coeff(:,2,1) = (K*(1.0_dp + K*nu_norm(:,1)**2.0_dp + 2.0_dp*(1.0_dp + K)*nu1(:,1)**4.0_dp - nu_norm(:,1)*(1.0_dp + K + (1.0_dp + 3.0_dp*K)*nu1(:,1)**2.0_dp) &
        !                 + nu1(:,1)**(2.0_dp)*(-1.0_dp + K + 2.0_dp*(1.0_dp + K)*nu2(:,1)**2.0_dp + 2.0_dp*(1.0_dp + K)*nu3(:,1)**2.0_dp)))/((1.0_dp + K)*(-1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        ! B_coeff(:,2,2) = (nu1(:,1)*(1.0_dp - 3.0_dp*K + (-1.0_dp + K)*K*nu_norm(:,1) + 2.0_dp*K*(-((-2.0_dp + K)*nu1(:,1)**2.0_dp) + nu2(:,1)**2.0_dp + nu3(:,1)**2.0_dp)))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! B_coeff(:,2,3) = (2.0_dp*K*(-1.0_dp + K*nu_norm(:,1) - (-1.0_dp + K)*nu1(:,1)**2.0_dp)*nu2(:,1))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! B_coeff(:,2,4) = (2.0_dp*K*(-1.0_dp + K*nu_norm(:,1) - (-1.0_dp + K)*nu1(:,1)**2.0_dp)*nu3(:,1))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! !### nu2(:,1) Equation Coefficients

        ! B_coeff(:,3,1) = (K*nu1(:,1)*nu2(:,1)*(-1.0_dp + K - (1.0_dp + 3.0_dp*K)*nu_norm(:,1) + 2.0_dp*(1.0_dp + K)*nu1(:,1)**2.0_dp + 2.0_dp*nu2(:,1)**2.0_dp &
        !                 + 2.0_dp*K*nu2(:,1)**2.0_dp + 2.0_dp*(1.0_dp + K)*nu3(:,1)**2.0_dp))/((1.0_dp + K)*(-1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        ! B_coeff(:,3,2) = (K*nu2(:,1)*(-1.0_dp - (1.0_dp + K)*nu_norm(:,1) + K*nu_norm(:,1)**2.0_dp - 2.0_dp*(-2.0_dp + K)*nu1(:,1)**2.0_dp &
        !                 + 2.0_dp*nu2(:,1)**2.0_dp + 2.0_dp*nu3(:,1)**2.0_dp))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! B_coeff(:,3,3) = nu1(:,1)*(-1.0_dp - (2.0_dp*(-1.0_dp + K)*K*nu2(:,1)**2.0_dp)/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp))*E11(:,1)

        ! B_coeff(:,3,4) = (-2.0_dp*(-1.0_dp + K)*K*nu1(:,1)*nu2(:,1)*nu3(:,1))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! !### nu3(:,1) Equation Coefficients

        ! B_coeff(:,4,1) = (K*nu1(:,1)*nu3(:,1)*(-1.0_dp + K - (1.0_dp + 3.0_dp*K)*nu_norm(:,1) + 2.0_dp*(1.0_dp + K)*nu1(:,1)**2.0_dp + 2.0_dp*nu2(:,1)**2.0_dp &
        !                 + 2.0_dp*K*nu2(:,1)**2.0_dp + 2.0_dp*(1.0_dp + K)*nu3(:,1)**2.0_dp))/((1.0_dp + K)*(-1.0_dp + K*nu_norm(:,1)))*E11(:,1)

        ! B_coeff(:,4,2) = (-(K*(1.0_dp + K*nu_norm(:,1)**2.0_dp - 2.0_dp*nu2(:,1)**2.0_dp + 2.0_dp*nu1(:,1)**2.0_dp*(-2.0_dp + K + 2.0_dp*K*(nu1(:,1)**(2.0_dp) + nu2(:,1)**2.0_dp)) &
        !                     + nu_norm(:,1)*(1.0_dp + K - 2.0_dp*K*(3.0_dp*nu1(:,1)**2.0_dp + nu2(:,1)**2.0_dp)))*nu3(:,1)) &
        !                     + 2.0_dp*K*(1.0_dp + K*nu_norm(:,1) - 2.0_dp*K*nu1(:,1)**2.0_dp)*nu3(:,1)**3.0_dp)/(-1.0_dp + K**(2.0_dp)*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! B_coeff(:,4,3) = (-2.0_dp*(-1.0_dp + K)*K*nu1(:,1)*nu2(:,1)*nu3(:,1))/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp)*E11(:,1)

        ! B_coeff(:,4,4) = nu1(:,1)*(-1.0_dp - (2.0_dp*(-1.0_dp + K)*K*nu3(:,1)**2.0_dp)/(-1.0_dp + K**2.0_dp*nu_norm(:,1)**2.0_dp))*E11(:,1)


    end subroutine compute_deriv_coefficients

    
    !########################################################################################################
    !# Compute RHS of PDE
    !########################################################################################################

    subroutine compute_rhs(t,u,dtu)
        implicit none

        ! Input/Output Variables
        real(dp), intent(in) :: t 
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar)
        ! real(dp), intent(in) :: K
        real(dp), intent(out) :: dtu(1-Ngz:Nx+Ngz,Nvar)
       

        ! System Variables:
        real(dp) :: Sigma_Minus(1-Ngz:Nx+Ngz,1), Sigma_Times(1-Ngz:Nx+Ngz,1), N_Minus(1-Ngz:Nx+Ngz,1), N_Times(1-Ngz:Nx+Ngz,1)
        real(dp) :: E11(1-Ngz:Nx+Ngz,1), Sigma2(1-Ngz:Nx+Ngz,1), Sigma3(1-Ngz:Nx+Ngz,1), LambdaTilde(1-Ngz:Nx+Ngz,1)
        real(dp) :: tau_log(1-Ngz:Nx+Ngz,1), tau(1-Ngz:Nx+Ngz,1)
        real(dp) :: mu(1-Ngz:Nx+Ngz,1), nu1(1-Ngz:Nx+Ngz,1), nu2(1-Ngz:Nx+Ngz,1), nu3(1-Ngz:Nx+Ngz,1)


        ! Auxiliary Variables: 
        real(dp) :: nu_norm(1-Ngz:Nx+Ngz,1)
        real(dp) :: Gamma2(1-Ngz:Nx+Ngz,1), T_00(1-Ngz:Nx+Ngz,1), T_01(1-Ngz:Nx+Ngz,1), T_02(1-Ngz:Nx+Ngz,1), T_03(1-Ngz:Nx+Ngz,1)
        real(dp) :: T_11(1-Ngz:Nx+Ngz,1), T_12(1-Ngz:Nx+Ngz,1), T_13(1-Ngz:Nx+Ngz,1), T_22(1-Ngz:Nx+Ngz,1), T_23(1-Ngz:Nx+Ngz,1)
        real(dp) :: T_33(1-Ngz:Nx+Ngz,1), Sigma_Plus(1-Ngz:Nx+Ngz,1), q(1-Ngz:Nx+Ngz,1), r(1-Ngz:Nx+Ngz,1)


        ! Flux/Reconstruction Variables:
        real(dp) :: cons_grav(1-Ngz:Nx+Ngz,4), prims_fluid(1-Ngz:Nx+Ngz,4)
        real(dp) :: cons_grav_plus(1-Ngz:Nx+Ngz,4), cons_grav_minus(1-Ngz:Nx+Ngz,4)
        real(dp) :: prims_fluid_plus(1-Ngz:Nx+Ngz,4), prims_fluid_minus(1-Ngz:Nx+Ngz,4)
        real(dp) :: E11_plus(1-Ngz:Nx+Ngz,1), E11_minus(1-Ngz:Nx+Ngz,1)
        real(dp) :: flux_grav(1-Ngz:Nx+Ngz,4), flux_grav_plus(1-Ngz:Nx+Ngz,4), flux_grav_minus(1-Ngz:Nx+Ngz,4)
        real(dp) :: B_fluid(1-Ngz:Nx+Ngz,4,4), B_fluid_plus(1-Ngz:Nx+Ngz,4,4),B_fluid_minus(1-Ngz:Nx+Ngz,4,4)
        real(dp) :: flux_grav_plus_E11(1-Ngz:Nx+Ngz,4), flux_grav_minus_E11(1-Ngz:Nx+Ngz,4)
        real(dp) :: char_speed(1-Ngz:Nx+Ngz,1)

        ! Non-Conservative Product Variables:
        real(dp) :: E11_deriv(1-Ngz:Nx+Ngz,1), E11_jump(1-Ngz:Nx+Ngz,1), tau_log_deriv(1-Ngz:Nx+Ngz,1), tau_log_jump(1-Ngz:Nx+Ngz,1)
        real(dp) :: nu1_deriv(1-Ngz:Nx+Ngz,1), nu1_jump(1-Ngz:Nx+Ngz,1), nu2_deriv(1-Ngz:Nx+Ngz,1), nu2_jump(1-Ngz:Nx+Ngz,1)
        real(dp) :: nu3_deriv(1-Ngz:Nx+Ngz,1), nu3_jump(1-Ngz:Nx+Ngz,1)
        real(dp) :: NC_deriv_grav(1-Ngz:Nx+Ngz,4), NC_jump_grav(1-Ngz:Nx+Ngz,4)
        real(dp) :: NC_deriv_fluid(1-Ngz:Nx+Ngz,4), NC_jump_fluid_plus(1-Ngz:Nx+Ngz,4), NC_jump_fluid_minus(1-Ngz:Nx+Ngz,4), NC_jump_fluid(1-Ngz:Nx+Ngz,4)

        ! Source Terms:
        real(dp) :: Sigma_Minus_source(1-Ngz:Nx+Ngz,1), Sigma_Times_source(1-Ngz:Nx+Ngz,1), N_Minus_source(1-Ngz:Nx+Ngz,1), N_Times_source(1-Ngz:Nx+Ngz,1)
        real(dp) :: tau_log_source(1-Ngz:Nx+Ngz,1), nu1_source(1-Ngz:Nx+Ngz,1), nu2_source(1-Ngz:Nx+Ngz,1), nu3_source(1-Ngz:Nx+Ngz,1)

        ! Time Derivatives:
        real(dp) :: dtSigma_Minus(1-Ngz:Nx+Ngz,1), dtSigma_Times(1-Ngz:Nx+Ngz,1), dtN_Minus(1-Ngz:Nx+Ngz,1), dtN_Times(1-Ngz:Nx+Ngz,1)
        real(dp) :: dtE11(1-Ngz:Nx+Ngz,1), dtSigma2(1-Ngz:Nx+Ngz,1), dtSigma3(1-Ngz:Nx+Ngz,1), dtLambdaTilde(1-Ngz:Nx+Ngz,1)
        real(dp) :: dttau_log(1-Ngz:Nx+Ngz,1), dtnu1(1-Ngz:Nx+Ngz,1), dtnu2(1-Ngz:Nx+Ngz,1), dtnu3(1-Ngz:Nx+Ngz,1)


        ! Re-scaling Parameter for Fluid Density
        real(dp) :: kappa

        !Temp variable
        real(dp) :: tmp_var(1-Ngz:Nx+Ngz,1)
        logical :: idx(1-Ngz:Nx+Ngz,1)


        ! Looping parameters
        integer :: i, j = 0

        ! ### Set re-scaling parameter
        ! kappa = (3.0_dp*K-1.0_dp)/(1.0_dp -K)
        kappa = 0.0_dp


        !# ---------------------------------
        !# Split the Solution Vector
        !# ---------------------------------

        Sigma_Minus(:,1) = u(:,1)
        Sigma_Times(:,1) = u(:,2)
        N_Minus(:,1) = u(:,3)
        N_Times(:,1) = u(:,4)
        E11(:,1) = u(:,5)
        Sigma2(:,1) = u(:,6)
        Sigma3(:,1) = u(:,7)
        LambdaTilde(:,1) = u(:,8)
        tau_log(:,1) = u(:,9)
        nu1(:,1) = u(:,10)
        nu2(:,1) = u(:,11)
        nu3(:,1) = u(:,12)

        ! Define Conserved Gravitational Vectors

        cons_grav(:,1) = Sigma_Minus(:,1)
        cons_grav(:,2) = Sigma_Times(:,1)
        cons_grav(:,3) = N_Minus(:,1)
        cons_grav(:,4) = N_Times(:,1)

        ! Define 'Primitive' Fluid Vector

        prims_fluid(:,1) = nu1(:,1)
        prims_fluid(:,2) = nu2(:,1)
        prims_fluid(:,3) = nu3(:,1)
        prims_fluid(:,4) = tau_log(:,1)


        ! In Vacuum, set fluid vars to zero:

        ! nu1(:,1) = 0.0_dp
        ! nu2(:,1) = 0.0_dp
        ! nu3(:,1) = 0.0_dp
        ! mu(:,1) = 0.0_dp


        !# --------------------------------------------
        !# Compute Auxiliary Variables
        !# --------------------------------------------

        !# Stress-Energy components are defined with *both* indices up!
        nu_norm = nu1**(2.0_dp) + nu2**(2.0_dp)  + nu3**(2.0_dp)
        Gamma2 = 1.0_dp/(1.0_dp - nu_norm)
        T_00 = exp(tau_log)
        T_01 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu1
        T_02 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu2
        T_03 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu3
        T_11 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu1*nu1 + K*T_00/(Gamma2*(1.0_dp+K*nu_norm))
        T_12 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu1*nu2 
        T_13 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu1*nu3 
        T_22 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu2*nu2 + K*T_00/(Gamma2*(1.0_dp+K*nu_norm))
        T_23 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu2*nu3 
        T_33 = (K+1)/(1.0_dp+K*nu_norm)*T_00*nu3*nu3 + K*T_00/(Gamma2*(1.0_dp+K*nu_norm))

        ! Sigma_Plus = 0.5_dp*(1.0_dp - Sigma_Times**(2.0_dp) - Sigma_Minus**(2.0_dp) - Sigma2**(2.0_dp) - Sigma3**(2.0_dp) &
        !             - N_Minus**(2.0_dp) - N_Times**(2.0_dp) - T_00 - LambdaTilde)

        ! q = 0.5_dp + (3.0_dp/2.0_dp)*(Sigma_Times**(2.0_dp) + Sigma_Minus**(2.0_dp) + Sigma2**(2.0_dp) + Sigma3**(2.0_dp)) - 3.0_dp*(Sigma2**(2.0_dp) + Sigma3**(2.0_dp)) &
        !     + (3.0_dp/2.0_dp)*(N_Minus**(2.0_dp) + N_Times**(2.0_dp)) + (3.0_dp/2.0_dp)*T_11 - (3.0_dp/2.0_dp)*LambdaTilde

        ! r = 3.0_dp*Sigma_Times*N_Minus - 3.0_dp*N_Times*Sigma_Minus - (3.0_dp/2.0_dp)*T_01


        !### With matter terms rescaled:
        Sigma_Plus = 0.5_dp*(1.0_dp - Sigma_Times**(2.0_dp) - Sigma_Minus**(2.0_dp) - Sigma2**(2.0_dp) - Sigma3**(2.0_dp) &
                    - N_Minus**(2.0_dp) - N_Times**(2.0_dp) - T_00*exp(kappa*t) - LambdaTilde)

        q = 0.5_dp + (3.0_dp/2.0_dp)*(Sigma_Times**(2.0_dp) + Sigma_Minus**(2.0_dp) + Sigma2**(2.0_dp) + Sigma3**(2.0_dp)) - 3.0_dp*(Sigma2**(2.0_dp) + Sigma3**(2.0_dp)) &
            + (3.0_dp/2.0_dp)*(N_Minus**(2.0_dp) + N_Times**(2.0_dp)) + (3.0_dp/2.0_dp)*T_11*exp(kappa*t) - (3.0_dp/2.0_dp)*LambdaTilde

        r = 3.0_dp*Sigma_Times*N_Minus - 3.0_dp*N_Times*Sigma_Minus - (3.0_dp/2.0_dp)*T_01*exp(kappa*t)


        !# --------------------------------------------
        !# Reconstruct at Cell Interfaces
        !# --------------------------------------------

        call linear_reconstruct(4,cons_grav,cons_grav_plus,cons_grav_minus)

        call linear_reconstruct(4,prims_fluid,prims_fluid_plus,prims_fluid_minus)

        call linear_reconstruct(1,E11,E11_plus,E11_minus)
            

        !# --------------------------------------------
        !# Compute the Characteristic Speeds
        !# --------------------------------------------

        call compute_characteristic_speeds(prims_fluid_plus,prims_fluid_minus,E11_plus,E11_minus,char_speed)


        !# --------------------------------------------
        !# Construct the Gravitational Fluxes
        !# --------------------------------------------

        ! Vacuum Char_Speed

        ! char_speed = max(abs(E11_minus),abs(E11_plus))

        ! CS_x = maxval(char_speed)

        call compute_flux_grav(cons_grav_plus,cons_grav_minus,flux_grav_plus,flux_grav_minus)
        

        ! Multiply fluxes by E11:
        do i = 1,4
            flux_grav_plus_E11(:,i) = E11_plus(:,1)*flux_grav_plus(:,i)
            flux_grav_minus_E11(:,i) = E11_minus(:,1)*flux_grav_minus(:,i)
        end do

        call local_lax_friedrichs(4,flux_grav_plus_E11,flux_grav_minus_E11,cons_grav_plus,cons_grav_minus,char_speed,flux_grav)



        !# --------------------------------------------------
        !# Construct the Fluid Derivative Coefficients
        !# --------------------------------------------------

        !# Cell-Average Coefficients

        call compute_deriv_coefficients(prims_fluid,E11,B_fluid)

        !# 'Plus' Reconstruction Coefficients

        call compute_deriv_coefficients(prims_fluid_plus,E11_plus,B_fluid_plus)

        !# 'Minus' Reconstruction Coefficients

        call compute_deriv_coefficients(prims_fluid_minus,E11_minus,B_fluid_minus)


        !# --------------------------------------------
        !# Construct the Non-Conservative Products
        !# --------------------------------------------


        !### First, we compute the gravitational non-conservative terms

        !### Compute the B(u)\del_{x}u terms

        E11_deriv(:,:) = 0.0_dp ! Initialise vector to zero

        ! E11_deriv at cell interfaces (size Nx)
        E11_deriv(1:Nx,1) = E11_minus(1:Nx,1) - E11_plus(0:Nx-1,1)

        !# Construct from reconstructed fluxes at cell interfaces 
        !# (i.e. taking an average of the fluxes at each side of an interface)
        NC_deriv_grav(:,:) = 0.0_dp


        NC_deriv_grav(1:Nx,1) = N_Times(1:Nx,1)*E11_deriv(1:Nx,1)
        NC_deriv_grav(1:Nx,2) = -N_Minus(1:Nx,1)*E11_deriv(1:Nx,1)
        NC_deriv_grav(1:Nx,3) = -Sigma_Times(1:Nx,1)*E11_deriv(1:Nx,1)
        NC_deriv_grav(1:Nx,4) = Sigma_Minus(1:Nx,1)*E11_deriv(1:Nx,1)
                

        !### Compute the non-conservative jump terms with linear segment path
        E11_jump(:,:) = 0.0_dp

        E11_jump = (E11_plus - E11_minus) !# This is the jump between reconstructions at the same interface

        do i = 1,4
            NC_jump_grav(:,i) = 0.5_dp*(0.5_dp*(flux_grav_plus(:,i)  + flux_grav_minus(:,i) )*E11_jump(:,1))
        end do



        !### Now, we compute the fluid non-conservative terms
        !### Importantly, since we have no flux terms in the fluid equations
        !### we *must* include the dissipation terms 
        !### 1/2*a_{i+1/2}*( w^{+}_{i+1/2} - w^{-}_{i+1/2} ) in the jump terms

        !### Compute the B(u)\del_{x}u terms

        tau_log_deriv(:,:) = 0.0_dp ! Initialise vector to zero
        nu1_deriv(:,:) = 0.0_dp 
        nu2_deriv(:,:) = 0.0_dp 
        nu3_deriv(:,:) = 0.0_dp 

        ! E11_deriv at cell interfaces (size Nx)
        nu1_deriv(1:Nx,1) = prims_fluid_minus(1:Nx,1) - prims_fluid_plus(0:Nx-1,1)
        nu2_deriv(1:Nx,1) = prims_fluid_minus(1:Nx,2) - prims_fluid_plus(0:Nx-1,2)
        nu3_deriv(1:Nx,1) = prims_fluid_minus(1:Nx,3) - prims_fluid_plus(0:Nx-1,3)
        tau_log_deriv(1:Nx,1) = prims_fluid_minus(1:Nx,4) - prims_fluid_plus(0:Nx-1,4)

        !# Construct from reconstructed fluxes at cell interfaces 
        !# (i.e. taking an average of the fluxes at each side of an interface)
        NC_deriv_fluid(:,:) = 0.0_dp


        NC_deriv_fluid(:,1) = B_fluid(:,1,1)*tau_log_deriv(:,1) + B_fluid(:,1,2)*nu1_deriv(:,1) + B_fluid(:,1,3)*nu2_deriv(:,1) + B_fluid(:,1,4)*nu3_deriv(:,1)
        NC_deriv_fluid(:,2) = B_fluid(:,2,1)*tau_log_deriv(:,1) + B_fluid(:,2,2)*nu1_deriv(:,1) + B_fluid(:,2,3)*nu2_deriv(:,1) + B_fluid(:,2,4)*nu3_deriv(:,1)
        NC_deriv_fluid(:,3) = B_fluid(:,3,1)*tau_log_deriv(:,1) + B_fluid(:,3,2)*nu1_deriv(:,1) + B_fluid(:,3,3)*nu2_deriv(:,1) + B_fluid(:,3,4)*nu3_deriv(:,1)
        NC_deriv_fluid(:,4) = B_fluid(:,4,1)*tau_log_deriv(:,1) + B_fluid(:,4,2)*nu1_deriv(:,1) + B_fluid(:,4,3)*nu2_deriv(:,1) + B_fluid(:,4,4)*nu3_deriv(:,1)
                

        !### Compute the non-conservative jump terms with linear segment path
        nu1_jump(:,:) = 0.0_dp
        nu2_jump(:,:) = 0.0_dp
        nu3_jump(:,:) = 0.0_dp
        tau_log_jump(:,:) = 0.0_dp

        nu1_jump(:,1) = (prims_fluid_plus(:,1) - prims_fluid_minus(:,1)) !# This is the jump between reconstructions at the same interface
        nu2_jump(:,1) = (prims_fluid_plus(:,2) - prims_fluid_minus(:,2))
        nu3_jump(:,1) = (prims_fluid_plus(:,3) - prims_fluid_minus(:,3))
        tau_log_jump(:,1) = (prims_fluid_plus(:,4) - prims_fluid_minus(:,4))

        !### log(tau) Jump Terms
        NC_jump_fluid(:,1) = 0.5_dp*( 0.5_dp*(B_fluid_plus(:,1,1)  + B_fluid_minus(:,1,1) )*tau_log_jump(:,1) &
                            +  0.5_dp*(B_fluid_plus(:,1,2)  + B_fluid_minus(:,1,2) )*nu1_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,1,3)  + B_fluid_minus(:,1,3) )*nu2_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,1,4)  + B_fluid_minus(:,1,4) )*nu3_jump(:,1))


        !### nu1 Jump Terms
        NC_jump_fluid(:,2) = 0.5_dp*( 0.5_dp*(B_fluid_plus(:,2,1)  + B_fluid_minus(:,2,1) )*tau_log_jump(:,1) &
                            +  0.5_dp*(B_fluid_plus(:,2,2)  + B_fluid_minus(:,2,2) )*nu1_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,2,3)  + B_fluid_minus(:,2,3) )*nu2_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,2,4)  + B_fluid_minus(:,2,4) )*nu3_jump(:,1)) 


        !### nu2 Jump Terms
        NC_jump_fluid(:,3) = 0.5_dp*( 0.5_dp*(B_fluid_plus(:,3,1)  + B_fluid_minus(:,3,1) )*tau_log_jump(:,1) &
                            +  0.5_dp*(B_fluid_plus(:,3,2)  + B_fluid_minus(:,3,2) )*nu1_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,3,3)  + B_fluid_minus(:,3,3) )*nu2_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,3,4)  + B_fluid_minus(:,3,4) )*nu3_jump(:,1))

        !### nu3 Jump Terms
        NC_jump_fluid(:,4) = 0.5_dp*( 0.5_dp*(B_fluid_plus(:,4,1)  + B_fluid_minus(:,4,1) )*tau_log_jump(:,1) &
                            +  0.5_dp*(B_fluid_plus(:,4,2)  + B_fluid_minus(:,4,2) )*nu1_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,4,3)  + B_fluid_minus(:,4,3) )*nu2_jump(:,1) & 
                            + 0.5_dp*(B_fluid_plus(:,4,4)  + B_fluid_minus(:,4,4) )*nu3_jump(:,1)) 
        



        !# --------------------------------------------
        !# Construct the Source Terms
        !# --------------------------------------------

        !### Gravitational Source Terms

        ! Sigma_Minus_source = (q + 3.0_dp*Sigma_Plus -2.0_dp)*Sigma_Minus + 2.0_dp*sqrt(3.0_dp)*Sigma_Times**(2.0_dp) + sqrt(3.0_dp)*(Sigma3**(2.0_dp) - Sigma2**(2.0_dp))&
        !                     - 2.0_dp*sqrt(3.0_dp)*N_Minus**(2.0_dp) + sqrt(3.0_dp)/2.0_dp*(T_22 - T_33)

        ! Sigma_Times_source = (q + 3.0_dp*Sigma_Plus -2.0_dp)*Sigma_Times - 2.0_dp*sqrt(3.0_dp)*N_Times*N_Minus - 2.0_dp*sqrt(3.0_dp)*Sigma_Times*Sigma_Minus &
        !                     + 2.0_dp*sqrt(3.0_dp)*Sigma2*Sigma3 + sqrt(3.0_dp)*T_23

        ! N_Minus_source = (q + 3.0_dp*Sigma_Plus)*N_Minus + 2.0_dp*sqrt(3.0_dp)*Sigma_Times*N_Times + 2.0_dp*sqrt(3.0_dp)*N_Minus*Sigma_Minus

        ! N_Times_source = (q + 3.0_dp*Sigma_Plus)*N_Times


        !### Gravitational Source terms with re-scaling of matter terms:
        Sigma_Minus_source = (q + 3.0_dp*Sigma_Plus -2.0_dp)*Sigma_Minus + 2.0_dp*sqrt(3.0_dp)*Sigma_Times**(2.0_dp) + sqrt(3.0_dp)*(Sigma3**(2.0_dp) - Sigma2**(2.0_dp))&
                            - 2.0_dp*sqrt(3.0_dp)*N_Minus**(2.0_dp) + sqrt(3.0_dp)/2.0_dp*(T_22 - T_33)*exp(kappa*t)

        Sigma_Times_source = (q + 3.0_dp*Sigma_Plus -2.0_dp)*Sigma_Times - 2.0_dp*sqrt(3.0_dp)*N_Times*N_Minus - 2.0_dp*sqrt(3.0_dp)*Sigma_Times*Sigma_Minus &
                            + 2.0_dp*sqrt(3.0_dp)*Sigma2*Sigma3 + sqrt(3.0_dp)*T_23*exp(kappa*t)

        N_Minus_source = (q + 3.0_dp*Sigma_Plus)*N_Minus + 2.0_dp*sqrt(3.0_dp)*Sigma_Times*N_Times + 2.0_dp*sqrt(3.0_dp)*N_Minus*Sigma_Minus

        N_Times_source = (q + 3.0_dp*Sigma_Plus)*N_Times





        !### Fluid Source Terms

        tau_log_source = (-1.0_dp - 3.0_dp*K + 2.0_dp*q + 2.0_dp*K*(1.0_dp + q)*nu_norm - nu3**(2.0_dp) - K*nu3**(2.0_dp) + sqrt(3.0_dp)*nu3**(2.0_dp)*Sigma_Minus + sqrt(3.0_dp)*K*nu3**(2.0_dp)*Sigma_Minus - &
                        (1.0_dp + K)*nu2**(2.0_dp)*(1.0_dp + sqrt(3.0_dp)*Sigma_Minus) + 3.0_dp*Sigma_Plus + 3.0_dp*K*Sigma_Plus + (1.0_dp + K)*nu1**(2.0_dp)*(-1.0_dp + 3.0_dp*Sigma_Plus) - &
                        2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*nu1*(nu3*Sigma2 + nu2*Sigma3) - 2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*nu2*nu3*Sigma_Times)/(1.0_dp + K*nu_norm)

        nu1_source = 2.0_dp*sqrt(3.0_dp)*N_Minus*nu2*nu3 + sqrt(3.0_dp)*N_Times*(-nu_norm + nu1**(2.0_dp) + 2.0_dp*nu3**(2.0_dp)) + &
                        (r - K*r - 2.0_dp*K**(2.0_dp)*r*nu_norm**(2.0_dp) + nu_norm*(K*(1.0_dp + K)*r + (1.0_dp + K)*(-1.0_dp + 3.0_dp*K)*nu1 + K*(1.0_dp + 5.0_dp*K)*r*nu1**(2.0_dp)) + &
                         nu1*(-((1.0_dp + K)*(-1.0_dp + 3.0_dp*K)) + r*nu1*(-1.0_dp + K - 2.0_dp*K**(2.0_dp) - 2.0_dp*K*(1.0_dp + K)*(nu1**(2.0_dp) + nu2**(2.0_dp) + nu3**(2.0_dp)))))/((1.0_dp + K)*(-1.0_dp + K*nu_norm)) + &
                        (sqrt(3.0_dp)*(-1.0_dp + K)*nu1*(nu_norm - nu1**(2.0_dp) - 2.0_dp*nu3**(2.0_dp))*Sigma_Minus)/(-1.0_dp + K*nu_norm) - &
                        (3.0_dp*(-1.0_dp + K)*nu1*(-1.0_dp + nu1**(2.0_dp))*Sigma_Plus)/(-1.0_dp + K*nu_norm) - &
                        (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K*nu_norm - (-1.0_dp + K)*nu1**(2.0_dp))*nu3*Sigma2)/(-1.0_dp + K*nu_norm) - &
                        (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K*nu_norm - (-1.0_dp + K)*nu1**(2.0_dp))*nu2*Sigma3)/(-1.0_dp + K*nu_norm) + &
                        (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu1*nu2*nu3*Sigma_Times)/(-1.0_dp + K*nu_norm)

        nu2_source = sqrt(3.0_dp)*N_Times*nu1*nu2 + (((1.0_dp + K)*(-1.0_dp + 3.0_dp*K)*(-1.0_dp + nu_norm) + r*(-1.0_dp + K - 2.0_dp*K**(2.0_dp) + K*(-1.0_dp + 3.0_dp*K)*nu_norm)*nu1)*nu2)/((1.0_dp + K)*(-1.0_dp + K*nu_norm)) + &
                    (sqrt(3.0_dp)*nu2*(1.0_dp - nu_norm - (-1.0_dp + K)*nu1**(2.0_dp) - 2.0_dp*(-1.0_dp + K)*nu3**(2.0_dp))*Sigma_Minus)/(-1.0_dp + K*nu_norm) + &
                    (3.0_dp*(K + nu1**(2.0_dp) - K*(nu_norm + nu1**(2.0_dp)))*nu2*Sigma_Plus)/(-1.0_dp + K*nu_norm) + &
                    (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu1*nu2*nu3*Sigma2)/(-1.0_dp + K*nu_norm) + (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu1*nu2**(2.0_dp)*Sigma3)/(-1.0_dp + K*nu_norm) + &
                    (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu2**(2.0_dp)*nu3*Sigma_Times)/(-1.0_dp + K*nu_norm)

        nu3_source = -2.0_dp*sqrt(3.0_dp)*N_Minus*nu1*nu2 - sqrt(3.0_dp)*N_Times*nu1*nu3 + &
                    (((1.0_dp+ K)*(-1.0_dp + 3.0_dp*K)*(-1.0_dp + nu_norm) + r*(-1.0_dp + K - 2.0_dp*K**(2.0_dp) + K*(-1.0_dp + 3.0_dp*K)*nu_norm)*nu1)*nu3)/((1.0_dp + K)*(-1.0_dp + K*nu_norm)) - &
                    (sqrt(3.0_dp)*nu3*(1.0_dp + (1.0_dp - 2.0_dp*K)*nu_norm + (-1.0_dp + K)*nu1**(2.0_dp) + 2.0_dp*(-1.0_dp + K)*nu3**(2.0_dp))*Sigma_Minus)/(-1.0_dp + K*nu_norm) + &
                    (3.0_dp*(K + nu1**(2.0_dp) - K*(nu_norm + nu1**(2.0_dp)))*nu3*Sigma_Plus)/(-1.0_dp + K*nu_norm) + (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu1*nu3**(2.0_dp)*Sigma2)/(-1.0_dp+ K*nu_norm) + &
                    (2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu1*nu2*nu3*Sigma3)/(-1.0_dp + K*nu_norm) - &
                    (2.0_dp*sqrt(3.0_dp)*nu2*(-1.0_dp + K*nu_norm - (-1.0_dp + K)*nu3**(2.0_dp))*Sigma_Times)/(-1.0_dp + K*nu_norm)



        ! tau_log_source = ((-1.0_dp - 3.0_dp*K + 2.0_dp*q + 2.0_dp*K*(1.0_dp + q)*nu_norm - nu3**2.0_dp - K*nu3**2.0_dp + sqrt(3.0_dp)*nu3**2.0_dp*Sigma_Minus + sqrt(3.0_dp)*K*nu3**2.0_dp*Sigma_Minus - &
        !                 (1.0_dp + K)*nu2**(2.0_dp)*(1.0_dp + sqrt(3.0_dp)*Sigma_Minus) + 3.0_dp*Sigma_Plus + 3.0_dp*K*Sigma_Plus + (1.0_dp + K)*nu1**(2.0_dp)*(-1.0_dp + 3.0_dp*Sigma_Plus) - &
        !                 2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*nu1*(nu3*Sigma2 + nu2*Sigma3) - 2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*nu2*nu3*Sigma_Times))/(1.0_dp + K*nu_norm)

        ! nu1_source = -(((-1.0_dp + K)*r + 2.0_dp*K**(2.0_dp)*r*nu_norm**(2.0_dp) + 2.0_dp*K*(1.0_dp + K)*r*nu1**4.0_dp + (1.0_dp + K)*nu1**3*(1.0_dp + K*(-3.0_dp + 4.0_dp*q) + (-3.0_dp + 9.0_dp*K)*Sigma_Plus) - &
        !             sqrt(3.0_dp)*(1.0_dp + K)*(-2.0_dp*N_Minus*nu2*nu3 + N_Times*(nu2**2.0_dp - nu3**2.0_dp) + 2.0_dp*nu3*Sigma2 + 2.0_dp*nu2*Sigma3) + &
        !             nu1**2.0_dp*(r + K*(-1.0_dp + 2.0_dp*K)*r + 2.0_dp*K*(1.0_dp + K)*r*nu3**2.0_dp - 2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K**2.0_dp)*nu3*Sigma2 + &
        !             2.0_dp*(1.0_dp + K)*nu2*(K*r*nu2 - sqrt(3.0_dp)*(-1.0_dp + K)*Sigma3)) + &
        !             K*nu_norm*(-((r + 5.0_dp*K*r)*nu1**2.0_dp) - 2.0_dp*(1.0_dp + K)*nu1*(2.0_dp*q + 3.0_dp*Sigma_Plus) + &
        !             (1.0_dp + K)*(-r + sqrt(3.0_dp)*(-2.0_dp*N_Minus*nu2*nu3 + N_Times*(nu2**2.0_dp - nu3**2.0_dp) + 2.0_dp*nu3*Sigma2 + 2.0_dp*nu2*Sigma3))) + &
        !             (1.0_dp + K)*nu1*(-1.0_dp + 3.0_dp*K - 3.0_dp*(-1.0_dp + K)*Sigma_Plus + nu2**2.0_dp*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q - sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus) + &
        !             nu3**2.0_dp*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q + sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus) - 2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu2*nu3*Sigma_Times))/&
        !             ((1.0_dp + K)*(-1.0_dp + K*nu_norm)))

        ! nu2_source = -((nu2*(2.0_dp*K*(1.0_dp + K)*r*nu1**3.0_dp + (1.0_dp + K)*nu1**2.0_dp*(1.0_dp + K*(-3.0_dp + 4.0_dp*q) + (-3.0_dp + 9.0_dp*K)*Sigma_Plus) - &
        !             K*nu_norm*((r + 5.0_dp*K*r + sqrt(3.0_dp)*(1.0_dp + K)*N_Times)*nu1 + (1.0_dp + K)*(4.0_dp*q - sqrt(3.0_dp)*Sigma_Minus + 3.0_dp*Sigma_Plus)) + &
        !             nu1*(r*(1.0_dp + K*(-1.0_dp + 2.0_dp*K) + 2.0_dp*K*(1.0_dp + K)*(nu2**2.0_dp + nu3**2.0_dp)) + &
        !             sqrt(3.0_dp)*(1.0_dp + K)*(N_Times - 2.0_dp*(-1.0_dp + K)*(nu3*Sigma2 + nu2*Sigma3))) + &
        !             (1.0_dp + K)*(-1.0_dp + 3.0_dp*K - sqrt(3.0_dp)*Sigma_Minus - 3.0_dp*K*Sigma_Plus + nu2**2.0_dp*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q - sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus) + &
        !             nu3**2.0_dp*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q + sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus) - 2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K)*nu2*nu3*Sigma_Times)))/&
        !             ((1.0_dp + K)*(-1.0_dp + K*nu_norm)))

        ! nu3_source = (-2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*N_Minus*(-1.0_dp + K*nu_norm)*nu1*nu2 - (1.0_dp + K)*nu3**3*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q + 2.0_dp*K*r*nu1 + sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus) + &
        !             nu3*(-2.0_dp*K*(1.0_dp + K)*r*nu1**3.0_dp - (1.0_dp + K)*nu1**2.0_dp*(1.0_dp + K*(-3.0_dp + 4.0_dp*q) + (-3.0_dp + 9.0_dp*K)*Sigma_Plus) + &
        !             (1.0_dp + K)*(1.0_dp - sqrt(3.0_dp)*Sigma_Minus + 3.0_dp*K*(-1.0_dp + Sigma_Plus) + K*nu_norm*(4.0_dp*q + sqrt(3.0_dp)*Sigma_Minus + 3.0_dp*Sigma_Plus) - &
        !             nu2**2.0_dp*(1.0_dp - 3.0_dp*K + 4.0_dp*K*q - sqrt(3.0_dp)*(-1.0_dp + K)*Sigma_Minus + 6.0_dp*K*Sigma_Plus)) + &
        !             nu1*(r*(-1.0_dp + K - 2.0_dp*K**2.0_dp + K*(1.0_dp + 5.0_dp*K)*nu_norm - 2.0_dp*K*(1.0_dp + K)*nu2**2.0_dp) - sqrt(3.0_dp)*(1.0_dp + K)*(N_Times*(-1.0_dp + K*nu_norm) - 2.0_dp*(-1.0_dp + K)*nu2*Sigma3))) - &
        !             2.0_dp*sqrt(3.0_dp)*(1.0_dp + K)*(-1.0_dp + K*nu_norm)*nu2*Sigma_Times + 2.0_dp*sqrt(3.0_dp)*(-1.0_dp + K**2.0_dp)*nu3**2.0_dp*(nu1*Sigma2 + nu2*Sigma_Times))/&
        !             ((1.0_dp + K)*(-1.0_dp + K*nu_norm))

        

        !# --------------------------------------------
        !# Update the PDE Variables
        !# --------------------------------------------
        
        !### Gravitational Equations
        dtSigma_Minus(:,:) = 0.0_dp
        dtSigma_Times(:,:) = 0.0_dp
        dtN_Minus(:,:) = 0.0_dp
        dtN_Times(:,:) = 0.0_dp
        

        dtSigma_Minus(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,1) - flux_grav(0:Nx-1,1) - NC_deriv_grav(1:Nx,1) &
                                - (NC_jump_grav(1:Nx,1) + NC_jump_grav(0:Nx-1,1))) &
                                + Sigma_Minus_source(1:Nx,1)
    
        dtSigma_Times(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,2) - flux_grav(0:Nx-1,2)  - NC_deriv_grav(1:Nx,2) &
                                    - (NC_jump_grav(1:Nx,2) + NC_jump_grav(0:Nx-1,2))) &
                                    + Sigma_Times_source(1:Nx,1)

        dtN_Minus(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,3) - flux_grav(0:Nx-1,3) - NC_deriv_grav(1:Nx,3) &
                                    - (NC_jump_grav(1:Nx,3)+ NC_jump_grav(0:Nx-1,3))) &
                                    + N_Minus_source(1:Nx,1)
        

        dtN_Times(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,4) - flux_grav(0:Nx-1,4) - NC_deriv_grav(1:Nx,4) &
                                    - (NC_jump_grav(1:Nx,4) + NC_jump_grav(0:Nx-1,4))) &
                                    + N_Times_source(1:Nx,1)



        !##### Finite Difference Version for Testing ##########

        ! tmp_var = FD_C4(N_Times)
        ! dtSigma_Minus(1:Nx,1) = -E11(1:Nx,1)*tmp_var(1:Nx,1) + Sigma_Minus_source(1:Nx,1)
        
        ! tmp_var = FD_C4(N_Minus)
        ! dtSigma_Times(1:Nx,1) = E11(1:Nx,1)*tmp_var(1:Nx,1)+ Sigma_Times_source(1:Nx,1)

        ! tmp_var = FD_C4(Sigma_Times)
        ! dtN_Minus(1:Nx,1) = E11(1:Nx,1)*tmp_var(1:Nx,1) + N_Minus_source(1:Nx,1)
        
        ! tmp_var = FD_C4(Sigma_Minus)
        ! dtN_Times(1:Nx,1) = -E11(1:Nx,1)*tmp_var(1:Nx,1) + N_Times_source(1:Nx,1)


        ! tmp_var = FD_C4(E11)
        ! dtSigma_Minus(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,1) - flux_grav(0:Nx-1,1)) + N_Times(1:Nx,1)*tmp_var(1:Nx,1) + Sigma_Minus_source(1:Nx,1)
        
        ! dtSigma_Times(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,2) - flux_grav(0:Nx-1,2)) - N_Minus(1:Nx,1)*tmp_var(1:Nx,1) + Sigma_Times_source(1:Nx,1)

        ! dtN_Minus(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,3) - flux_grav(0:Nx-1,3)) - Sigma_Times(1:Nx,1)*tmp_var(1:Nx,1)+ N_Minus_source(1:Nx,1)
        
        ! dtN_Times(1:Nx,1) = -1.0_dp/(dx)*(flux_grav(1:Nx,4) - flux_grav(0:Nx-1,4)) + Sigma_Minus(1:Nx,1)*tmp_var(1:Nx,1) + N_Times_source(1:Nx,1)


        !### Fluid Equations
        dttau_log(:,:) = 0.0_dp
        dtnu1(:,:) = 0.0_dp
        dtnu2(:,:) = 0.0_dp
        dtnu3(:,:) = 0.0_dp


        ! dttau_log(1:Nx,1) = -1.0_dp/(dx)*(- NC_deriv_fluid(1:Nx,1) - (NC_jump_fluid_minus(1:Nx,1) + NC_jump_fluid_plus(0:Nx-1,1))) &
        !                         + tau_log_source(1:Nx,1)
    
        ! dtnu1(1:Nx,1) = -1.0_dp/(dx)*(- NC_deriv_fluid(1:Nx,2) - (NC_jump_fluid_minus(1:Nx,2)+ NC_jump_fluid_plus(0:Nx-1,2))) &
        !                             + nu1_source(1:Nx,1)

        ! dtnu2(1:Nx,1) = -1.0_dp/(dx)*(- NC_deriv_fluid(1:Nx,3) - (NC_jump_fluid_minus(1:Nx,3) + NC_jump_fluid_plus(0:Nx-1,3))) &
        !                             + nu2_source(1:Nx,1)
        

        ! dtnu3(1:Nx,1) = -1.0_dp/(dx)*(- NC_deriv_fluid(1:Nx,4) - (NC_jump_fluid_minus(1:Nx,4) + NC_jump_fluid_plus(0:Nx-1,4))) &
        !                             + nu3_source(1:Nx,1)



        dttau_log(1:Nx,1) = -1.0_dp/(dx)*(0.5_dp*char_speed(1:Nx,1)*tau_log_jump(1:Nx,1) - 0.5_dp*char_speed(0:Nx-1,1)*tau_log_jump(0:Nx-1,1) &
                                - NC_deriv_fluid(1:Nx,1) - (NC_jump_fluid(1:Nx,1) + NC_jump_fluid(0:Nx-1,1))) &
                                + tau_log_source(1:Nx,1)
    
        dtnu1(1:Nx,1) = -1.0_dp/(dx)*(0.5_dp*char_speed(1:Nx,1)*nu1_jump(1:Nx,1) - 0.5_dp*char_speed(0:Nx-1,1)*nu1_jump(0:Nx-1,1) &
                                - NC_deriv_fluid(1:Nx,2) - (NC_jump_fluid(1:Nx,2) + NC_jump_fluid(0:Nx-1,2))) &
                                + nu1_source(1:Nx,1)

        dtnu2(1:Nx,1) = -1.0_dp/(dx)*(0.5_dp*char_speed(1:Nx,1)*nu2_jump(1:Nx,1) - 0.5_dp*char_speed(0:Nx-1,1)*nu2_jump(0:Nx-1,1) &
                                - NC_deriv_fluid(1:Nx,3) - (NC_jump_fluid(1:Nx,3) + NC_jump_fluid(0:Nx-1,3))) &
                                + nu2_source(1:Nx,1)
        

        dtnu3(1:Nx,1) = -1.0_dp/(dx)*(0.5_dp*char_speed(1:Nx,1)*nu3_jump(1:Nx,1) - 0.5_dp*char_speed(0:Nx-1,1)*nu3_jump(0:Nx-1,1) &
                                - NC_deriv_fluid(1:Nx,4) - (NC_jump_fluid(1:Nx,4) + NC_jump_fluid(0:Nx-1,4))) &
                                + nu3_source(1:Nx,1)



        !# For future evolution swap sign of dissipation terms:

        ! dttau_log(1:Nx,1) = -1.0_dp/(dx)*(-0.5_dp*char_speed(1:Nx,1)*tau_log_jump(1:Nx,1) + 0.5_dp*char_speed(0:Nx-1,1)*tau_log_jump(0:Nx-1,1) &
        !                                 - NC_deriv_fluid(1:Nx,1) - (NC_jump_fluid(1:Nx,1) + NC_jump_fluid(0:Nx-1,1))) &
        !                                 + tau_log_source(1:Nx,1)

        ! dtnu1(1:Nx,1) = -1.0_dp/(dx)*(-0.5_dp*char_speed(1:Nx,1)*nu1_jump(1:Nx,1) + 0.5_dp*char_speed(0:Nx-1,1)*nu1_jump(0:Nx-1,1) &
        !                         - NC_deriv_fluid(1:Nx,2) - (NC_jump_fluid(1:Nx,2) + NC_jump_fluid(0:Nx-1,2))) &
        !                         + nu1_source(1:Nx,1)

        ! dtnu2(1:Nx,1) = -1.0_dp/(dx)*(-0.5_dp*char_speed(1:Nx,1)*nu2_jump(1:Nx,1) + 0.5_dp*char_speed(0:Nx-1,1)*nu2_jump(0:Nx-1,1) &
        !                         - NC_deriv_fluid(1:Nx,3) - (NC_jump_fluid(1:Nx,3) + NC_jump_fluid(0:Nx-1,3))) &
        !                         + nu2_source(1:Nx,1)
        

        ! dtnu3(1:Nx,1) = -1.0_dp/(dx)*(-0.5_dp*char_speed(1:Nx,1)*nu3_jump(1:Nx,1) + 0.5_dp*char_speed(0:Nx-1,1)*nu3_jump(0:Nx-1,1) &
        !                         - NC_deriv_fluid(1:Nx,4) - (NC_jump_fluid(1:Nx,4) + NC_jump_fluid(0:Nx-1,4))) &
        !                         + nu3_source(1:Nx,1)



        !####### CFL type condition for improved positivity preservation #########

        ! ODE Version:
        ! CS_source = minval((1.0_dp - 1.0e-14_dp - sqrt(nu_norm(1:Nx,1)))/(sqrt(nu1_source(1:Nx,1)**(2.0_dp) + nu2_source(1:Nx,1)**(2.0_dp)  + nu3_source(1:Nx,1)**(2.0_dp) )))

        ! PDE Version:
        ! CS_source = minval((1.0_dp - 0.0e-14_dp - sqrt(nu_norm(1:Nx,1)))/(sqrt(dtnu1(1:Nx,1)**(2.0_dp) + dtnu2(1:Nx,1)**(2.0_dp)  + dtnu3(1:Nx,1)**(2.0_dp) )))


        ! print *, nu_norm

        !# --------------------------------------------
        !# Update the ODE Variables
        !# --------------------------------------------
        dtE11(:,:) = 0.0_dp
        dtSigma2(:,:) = 0.0_dp
        dtSigma3(:,:) = 0.0_dp
        dtLambdaTilde(:,:) = 0.0_dp

        ! dtE11 = (q + 3.0_dp*Sigma_Plus)*E11

        ! dtSigma2 = (q - 2.0_dp + sqrt(3.0_dp)*Sigma_Minus)*Sigma2 - 2.0_dp*sqrt(3.0_dp)*Sigma3*Sigma_Times + sqrt(3.0_dp)*T_13

        ! dtSigma3 = (q - 2.0_dp - sqrt(3.0_dp)*Sigma_Minus)*Sigma3 + sqrt(3.0_dp)*T_12

        ! dtLambdaTilde = 2.0_dp*(q + 1.0_dp)*LambdaTilde


        !### ODEs with Matter terms rescaled:
        dtE11 = (q + 3.0_dp*Sigma_Plus)*E11

        dtSigma2 = (q - 2.0_dp + sqrt(3.0_dp)*Sigma_Minus)*Sigma2 - 2.0_dp*sqrt(3.0_dp)*Sigma3*Sigma_Times + sqrt(3.0_dp)*T_13*exp(kappa*t)

        dtSigma3 = (q - 2.0_dp - sqrt(3.0_dp)*Sigma_Minus)*Sigma3 + sqrt(3.0_dp)*T_12*exp(kappa*t)

        dtLambdaTilde = 2.0_dp*(q + 1.0_dp)*LambdaTilde



        ! if (j==180) then
        !     print *, t
        !     open(unit=10, file="test.dat", status="replace", action="write", &
        !     access="stream", form="unformatted")
        !     write(10) dtSigma_Times
        !     flush(10)
        !     close(10)
        ! end if

        ! j = j+1

        ! print *, j


        !# --------------------------------------------
        !# Collect Time Derivatives
        !# --------------------------------------------

        dtu(:,1) = dtSigma_Minus(:,1)
        dtu(:,2) = dtSigma_Times(:,1)
        dtu(:,3) = dtN_Minus(:,1)
        dtu(:,4) = dtN_Times(:,1)
        dtu(:,5) = dtE11(:,1)
        dtu(:,6) = dtSigma2(:,1)
        dtu(:,7) = dtSigma3(:,1)
        dtu(:,8) = dtLambdaTilde(:,1)
        dtu(:,9) = dttau_log(:,1)
        dtu(:,10) = dtnu1(:,1)
        dtu(:,11) = dtnu2(:,1)
        dtu(:,12) = dtnu3(:,1)
        

    end subroutine compute_rhs

end module T2_fluid_mod