program evolve_system
    use real_type_mod, only: dp
    use time_stepping_mod, only: rk3, rk4, rk1
    use boundary_conditions_mod, only: periodic_bc
    use T2_fluid_mod, only: CS_source, CS_x
    use read_parameters_mod
    use finite_difference_mod
    use hdf5_IO_mod
    implicit none

    integer :: beginning, steps, end
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: CFL, dt, rate, tsave
    

    !# Initial Data Variables
    real(dp), allocatable :: y0(:,:)
    real(dp), allocatable :: Sigma_Minus(:,:), Sigma_Times(:,:), N_Minus(:,:), N_Times(:,:) 
    real(dp), allocatable :: E11(:,:), Sigma2(:,:), Sigma3(:,:), LambdaTilde(:,:)
    real(dp), allocatable :: tau(:,:), S1(:,:), S2(:,:), S3(:,:), tau_log(:,:)
    real(dp), allocatable :: nu1(:,:), nu2(:,:), nu3(:,:), mu(:,:), Gamma2(:,:)
    real(dp), allocatable :: r(:,:), CM1(:,:), CM2(:,:), CBeta(:,:)
    real(dp) :: a_param, b_param, c_param, d_param


    !# File and Storage Buffer Variables
    integer, parameter :: buffer_size = 1000
    real(dp), allocatable :: solution_buffer(:,:,:)
    real(dp) :: time_buffer(buffer_size)
    integer :: buffer_index = 0, ts=0

    call system_clock(beginning, rate)


    !#######################################
    !# Load Simulation Parameters
    !#######################################
    call read_params()

    ! Initial timestep size
    CFL = 0.01_dp
    ! CFL = 1.0e-4_dp
    dt = CFL*dx

    ! Allocate I/O Buffer
    allocate(solution_buffer(Nx, Nvar, buffer_size))


    ! Initialise HDF5 File
    call initialise_hdf(trim(filename),Nx,500,coordinates_x,K)


    !#########################################
    ! Set Initial Data
    !########################################

    !# Allocate Initial Data Arrays
    allocate(y0(1-Ngz:Nx+Ngz,Nvar))
    allocate(nu1(1-Ngz:Nx+Ngz,1))
    allocate(nu2(1-Ngz:Nx+Ngz,1))
    allocate(nu3(1-Ngz:Nx+Ngz,1))
    allocate(Gamma2(1-Ngz:Nx+Ngz,1))
    allocate(mu(1-Ngz:Nx+Ngz,1))
    allocate(tau(1-Ngz:Nx+Ngz,1))
    allocate(tau_log(1-Ngz:Nx+Ngz,1))
    allocate(S1(1-Ngz:Nx+Ngz,1))
    allocate(S2(1-Ngz:Nx+Ngz,1))
    allocate(S3(1-Ngz:Nx+Ngz,1))
    allocate(Sigma_Minus(1-Ngz:Nx+Ngz,1))
    allocate(Sigma_Times(1-Ngz:Nx+Ngz,1))
    allocate(N_Minus(1-Ngz:Nx+Ngz,1))
    allocate(N_Times(1-Ngz:Nx+Ngz,1))
    allocate(E11(1-Ngz:Nx+Ngz,1))
    allocate(Sigma2(1-Ngz:Nx+Ngz,1))
    allocate(Sigma3(1-Ngz:Nx+Ngz,1))
    allocate(LambdaTilde(1-Ngz:Nx+Ngz,1))
    
    ! -------------------------------------
    ! `Generic' Cosmology Initial Data
    ! -------------------------------------

    !### Define perturbation parameters
    a_param = -1.0_dp
    b_param = -1.0_dp
    c_param = 0.01_dp
    d_param = 0.1_dp


    !####### Fluid Variables
    nu1(:,1) = c_param*sin(coordinates_x)

    nu2(:,1) = (-3.0_dp*a_param*c_param)/(2.0_dp*sqrt(3.0_dp))*sin(coordinates_x)

    nu3(:,1) = (-3.0_dp*b_param*c_param)/(2.0_dp*sqrt(3.0_dp))*sin(coordinates_x)

    mu(:,1) = 0.0_dp*coordinates_x + 1.0_dp

    Gamma2(:,1) = 1.0_dp/(1.0_dp - nu1(:,1)*nu1(:,1) - nu2(:,1)*nu2(:,1) - nu3(:,1)*nu3(:,1))

    !####### Define Conserved Variables

    tau(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1) - K*mu(:,1)
    S1(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu1(:,1)
    S2(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu2(:,1)
    S3(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu3(:,1)

    !### Log T^{00} variable:
    tau_log = log(tau)

    !####### Cosmological Constant Variable

    LambdaTilde(:,1) = 0.0_dp*coordinates_x + 0.0_dp

    !####### Gravitational Variables

    E11(:,1) = 0.0_dp*coordinates_x + 1.0_dp

    Sigma_Minus(:,1) = d_param*sin(coordinates_x) 

    Sigma_Times(:,1) = d_param*cos(coordinates_x)

    N_Minus(:,1) = 0.0_dp*coordinates_x 

    N_Times(:,1) = 0.0_dp*coordinates_x

    Sigma3(:,1) = 0.0_dp*coordinates_x + a_param

    Sigma2(:,1) = 0.0_dp*coordinates_x + b_param

    ! -------------------------------------
    ! Future Asymptotics Test Initial Data
    ! -------------------------------------

    !####### Fluid Variables
    ! nu1(:,1) = 0.1_dp*sin(coordinates_x)
    ! nu2(:,1) = 0.0_dp*sin(coordinates_x)
    ! nu3(:,1) = 0.0_dp*sin(coordinates_x)
    ! mu(:,1) = 0.0_dp*coordinates_x + 0.1_dp

    ! Gamma2(:,1) = 1.0_dp/(1.0_dp - nu1(:,1)*nu1(:,1) - nu2(:,1)*nu2(:,1) - nu3(:,1)*nu3(:,1))

    ! !####### Define Conserved Variables

    ! tau(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1) - K*mu(:,1)
    ! S1(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu1(:,1)
    ! S2(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu2(:,1)
    ! S3(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu3(:,1)

    ! !### Log T^{00} variable:
    ! tau_log = log(tau)

    ! !####### Cosmological Constant Variables

    ! LambdaTilde(:,1) = 0.0_dp*coordinates_x + 0.1_dp

    ! !####### Gravitational Variables

    ! E11(:,1) = 0.0_dp*coordinates_x + 1.0_dp

    ! Sigma_Minus(:,1) = 0.0_dp*coordinates_x 

    ! Sigma_Times(:,1) = S1(:,1) 

    ! N_Minus(:,1) = 0.0_dp*coordinates_x + 0.5_dp

    ! N_Times(:,1) = 0.0_dp*coordinates_x

    ! Sigma3(:,1) = 0.0_dp*coordinates_x

    ! Sigma2(:,1) = 0.0_dp*coordinates_x

    ! -------------------------------------
    ! Exact Vacuum Gowdy Test Solution
    ! -------------------------------------

    ! !####### Fluid Variables
    ! nu1(:,1) = 0.0_dp*sin(coordinates_x)
    ! nu2(:,1) = 0.0_dp*sin(coordinates_x)
    ! nu3(:,1) = 0.0_dp*sin(coordinates_x)
    ! mu(:,1) = 0.0_dp*coordinates_x + 0.0_dp

    ! Gamma2(:,1) = 1.0_dp/(1.0_dp - nu1(:,1)*nu1(:,1) - nu2(:,1)*nu2(:,1) - nu3(:,1)*nu3(:,1))

    ! ! !####### Define Conserved Variables

    ! tau(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1) - K*mu(:,1)
    ! S1(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu1(:,1)
    ! S2(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu2(:,1)
    ! S3(:,1) = (K+1.0_dp)*Gamma2(:,1)*mu(:,1)*nu3(:,1)

    ! ! !####### Cosmological Constant Variables

    ! LambdaTilde(:,1) = 0.0_dp*coordinates_x + 0.0_dp

    ! ! !####### Gravitational Variables

    ! E11(:,1) = 0.0_dp*coordinates_x + 1.0_dp

    ! Sigma_Minus(:,1) = 0.0_dp*coordinates_x - 1/(2.0_dp*sqrt(3.0_dp))

    ! Sigma_Times(:,1) = 0.0_dp*coordinates_x + 1/(sqrt(3.0_dp))*exp(t0)*sin(exp(2.0_dp*t0)-2.0_dp*coordinates_x)

    ! N_Minus(:,1) = 0.0_dp*coordinates_x - 1/(sqrt(3.0_dp))*exp(t0)*sin(exp(2.0_dp*t0)-2.0_dp*coordinates_x)

    ! N_Times(:,1) = 0.0_dp*coordinates_x

    ! Sigma3(:,1) = 0.0_dp*coordinates_x

    ! Sigma2(:,1) = 0.0_dp*coordinates_x


    !###################### Check Initial Constraint Violation ###################
    ! allocate(r(1-Ngz:Nx+Ngz,1))
    ! allocate(CM1(1-Ngz:Nx+Ngz,1))
    ! allocate(CM2(1-Ngz:Nx+Ngz,1))
    ! allocate(CBeta(1-Ngz:Nx+Ngz,1))


    ! r = 3.0_dp*Sigma_Times*N_Minus - 3.0_dp*N_Times*Sigma_Minus - 3.0_dp/2.0_dp*S1

    ! CM1 = E11*FD_C4(Sigma3) - r*Sigma3 - sqrt(3.0_dp)*Sigma3*N_Times + sqrt(3.0_dp)*S2

    ! CM2 = E11*FD_C4(Sigma2) - r*Sigma2 - sqrt(3.0_dp)*Sigma2*N_Times + sqrt(3.0_dp)*S3 + 2.0_dp*sqrt(3.0_dp)*Sigma3*N_Minus

    ! CBeta = E11*FD_C4(LambdaTilde) - 2*LambdaTilde*r


    ! print *, CM1



    !# Collect Initial Data
    y0(:,1) = Sigma_Minus(:,1)
    y0(:,2) = Sigma_Times(:,1)
    y0(:,3) = N_Minus(:,1)
    y0(:,4) = N_Times(:,1)
    y0(:,5) = E11(:,1)
    y0(:,6) = Sigma2(:,1)
    y0(:,7) = Sigma3(:,1)
    y0(:,8) = LambdaTilde(:,1)
    y0(:,9) = tau_log(:,1)
    y0(:,10) = nu1(:,1)
    y0(:,11) = nu2(:,1)
    y0(:,12) = nu3(:,1)

    !# Apply Boundary Conditions to Initial Data
    call periodic_bc(Nx,Ngz,Nvar,y0)

    !#########################################
    ! Set up files for saving data
    !########################################

    ! Save initial data to buffer arrays
    buffer_index = 1
    solution_buffer(:, :, buffer_index) = y0(1:Nx,:)
    time_buffer(buffer_index) = t0

    
    CS_x = 1.0_dp
    CS_source = 1.0_dp
    !#########################################
    ! Evolution Routine
    !########################################
    steps = 0
    tsave = -0.01_dp ! Past
    ! tsave = 0.01_dp ! Future
    ! do while (t0<tend) ! Future Evolution
    do while (t0>tend)  ! Past Evolution

        steps = steps + 1

        ! Adjust timestep
        ! dt = min(CFL*dx/CS_x,CFL*dx) ! To the future
        dt = -min(CFL*dx/CS_x,CFL*dx) ! To the past
        ! dt = -min(CFL*dx/CS_x,CFL*dx,1.0e-4_dp) 
        ! dt = -min(CFL*dx/CS_x, CFL*dx, CS_source*CFL) ! To the past

        call rk4(t0,y0,dt)
        ! call rk1_imp(t0,y0,dt)
        t0 = t0 + dt


        ! print *, "t is", t0
        ! print *, "dt is", dt
        ! print *, "CS_x is", CS_x
        ! print *, "CS_source is", CS_source
        ! print *, " "

        ! exit

        ! if (abs(dt) <= 1.0e-15_dp) then
        !     print *, "Here!"
        !     exit
        ! end if


        !# Save solution data every 100 timesteps
        ! if (steps == 5) then 
        !     buffer_index = buffer_index + 1
        !     solution_buffer(:, :, buffer_index) = y0(1:Nx, :)
        !     time_buffer(buffer_index) = t0
        !     ! print *, t0

        !     ! if (buffer_index == 60 .and. ts == 0) then
        !     !     ts = 1
        !     !     open(unit=10, file="test.dat", status="replace", action="write", &
        !     !     access="stream", form="unformatted")
        !     !     write(10) y0(:, 2)
        !     !     flush(10)
        !     !     close(10)
        !     ! end if
            
        !     ! Save buffer array to file when full
        !     if (buffer_index == buffer_size) then
        !         call write_to_hdf(solution_buffer,time_buffer,buffer_index,Nx)
        !         buffer_index = 0
        !     end if
            
        !     steps = 0
        ! end if

        !# Past Saving
        if (t0 <= tsave) then 
            buffer_index = buffer_index + 1
            solution_buffer(:, :, buffer_index) = y0(1:Nx, :)
            time_buffer(buffer_index) = t0
            tsave = tsave - 0.01_dp
            print *, "t is", t0
            print *, "dt is", dt
            print *, " "

            ! if (buffer_index == 60 .and. ts == 0) then
            !     ts = 1
            !     open(unit=10, file="test.dat", status="replace", action="write", &
            !     access="stream", form="unformatted")
            !     write(10) y0(:, 2)
            !     flush(10)
            !     close(10)
            ! end if
            
            ! Save buffer array to file when full
            if (buffer_index == buffer_size) then
                call write_to_hdf(solution_buffer,time_buffer,buffer_index,Nx)
                buffer_index = 0
            end if
            
            steps = 0
        end if

        !# Future Saving
        ! if (t0 >= tsave) then 
        !     buffer_index = buffer_index + 1
        !     solution_buffer(:, :, buffer_index) = y0(1:Nx, :)
        !     time_buffer(buffer_index) = t0
        !     tsave = tsave + 0.01_dp
        !     ! print *, "t is", t0
        !     ! print *, "dt is", dt
        !     ! print *, " "

        !     ! if (buffer_index == 60 .and. ts == 0) then
        !     !     ts = 1
        !     !     open(unit=10, file="test.dat", status="replace", action="write", &
        !     !     access="stream", form="unformatted")
        !     !     write(10) y0(:, 2)
        !     !     flush(10)
        !     !     close(10)
        !     ! end if
            
        !     ! Save buffer array to file when full
        !     if (buffer_index == buffer_size) then
        !         call write_to_hdf(solution_buffer,time_buffer,buffer_index,Nx)
        !         buffer_index = 0
        !     end if
            
        !     steps = 0
        ! end if

    end do

    ! Save any remaining data in buffer array
    if (buffer_index > 0) then
        call write_to_hdf(solution_buffer,time_buffer,buffer_index,Nx)
        call close_hdf()
    end if

    call system_clock(end)
    print *, "elapsed time: ", real(end - beginning) / real(rate)



contains 
    

end program evolve_system