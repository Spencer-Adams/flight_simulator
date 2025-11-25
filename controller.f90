module sim_m
    use adams_m
    use controller_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m
    implicit none
    ! Variables within sim_m
    real :: mass
    real :: I(3,3)
    real :: Ixxb, Iyyb, Izzb, Ixyb, Ixzb, Iyzb
    real :: Iinv(3,3)
    real :: h_gyro(3,3)
    real :: hdot_gyro(3)
    real :: hx, hy, hz ! non-zero terms of h_gyro
    real :: FM(6)
    real :: y_init(13)
    real :: controls(4)
    real :: controls_from_connect(4)
    real :: V_initial, alpha_initial, beta_initial
    real :: rho0,Z_temp,T_temp,P_temp,a_temp,mu_temp
    real :: trim_state(13)
    logical :: rk4_verbose, is_trim_sideslip_angle, is_use_controls
    real :: trim_elevation_angle, trim_sideslip_angle, trim_bank_angle
    real :: trim_azimuth_angle, p_wind, trim_climb_angle

    type(connection) :: my_controller 
    type(json_value), pointer :: j_main

    ! aero coefficients 
    real, allocatable :: aero_ref_location(:), eul0(:), angular_rates(:) 
    character(len=:),allocatable :: init_type, trim_type
    character(len=:),allocatable :: is_elevation_or_climb, is_bank_or_beta_for_shss
    real :: finite_diff_step, relax_factor, newton_tol
    real :: trim_array(9)
    real :: sref, long_ref, lat_ref
    real :: CL0, CLa, CLahat, CLqbar, CLde
    real :: CDL0, CDL1, CDL2, CDS2, CDqbar, CDaqbar, CDde, CDade, CDde2
    real :: CSb, CSpbar, CSapbar, CSrbar, CSda, CSdr
    real :: Clb, Clpbar, Clrbar, Clarbar, Clda, Cldr
    real :: Cm0, Cma, Cmqbar, Cmahat, Cmde 
    real :: Cnb, Cnpbar, Cnapbar, Cnrbar, Cnda, Cnada, Cndr 
    real :: Thrust0, Ta
    real :: weight !! might switch this out for more general stuff.
    real ::  t, u_0, v_0, w_0, p_0, q_0, r_0, x_0, y_0, z_0
    real :: dt, tf, delta_t_over_2, delta_t_over_6
    integer :: n, newton_max_iter
    contains
    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        real :: y(13), y1(13), s(14)
        ! real :: Z_temp, P_temp, T_temp, a_temp, mu_temp
        real :: cpu_start_time, cpu_end_time, time1, time2, actual_time, integrated_time
        integer :: io_unit
        logical :: real_time
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        n = 13
        ! call std_atm_English(0.0,Z_temp,T_temp,P_temp,rho0,a_temp,mu_temp)
        ! initial conditions 
        t = 0.0
        y = y_init
        real_time = .false.
        if (abs(dt)<TOLERANCE) then 
            real_time = .true.
            time1 = get_time()
            y1 = runge_kutta(t, y,dt)
            call quat_norm(y1(10:13))
            time2 = get_time()
            dt = time2-time1 
            write(*,*) "Dt = ", dt
            if (dt == 0.0) then 
                write(*,*) "THE DT IS EXACTLY ZERO!!!"
            end if 
            y = y_init 
            t = 0.0
        end if
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        open(newunit=io_unit, file='output.txt', status='replace', action='write')
        call quat_norm(y(10:13))
        write(io_unit,*) "        t[s]                   u[ft/s]               v[ft/s]"// &
                   "                 w[ft/s]                p[rad/s]             q[rad/s]"  // &
                   "            r[rad/s]                x[ft]                  y[ft]" // &                   
                   "                  z[ft]                    e0                   ex" // &
                   "                      ey                   ez"
        write(io_unit,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        ! write(*,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        ! do while(y(9)<0.0 .and. t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
        cpu_start_time = get_time()
        time1 = cpu_start_time
        integrated_time = 0.0

        do while(t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
            if (is_use_controls) then 
                controls = connect_controls%recv()
            end if
            y1 = runge_kutta(t,y,dt)
            call quat_norm(y1(10:13))
            y = y1
            t = t + dt 
            integrated_time = integrated_time + dt
            ! write(io_unit,'(14E22.13)') t,y(:)
            s(1) = t
            s(2:14) = y(1:13)
            ! write(*,'(14E22.13)') s
            call graphics%send(s)
            ! write(*,'(14E22.13)') t,y(:)
            if(real_time) then ! get estimate for next dt time step
                time2 = get_time()
                dt = time2-time1 
                time1 = time2 
            end if 

        end do
        cpu_end_time = get_time()
        actual_time = cpu_end_time-cpu_start_time
        write(*,*) '    Total integrated time [s] = ', integrated_time
        write(*,*) 'Total actual elapsed time [s]= ', actual_time
        write(*,*) 'Total error in time [s] = ', integrated_time - actual_time
    end subroutine run


    subroutine init(filename)
        implicit none 
        character(100), intent(in) :: filename
        type(json_value), pointer :: j_connections, j_graphics, j_controls
        ! type2, intent(out) ::  arg2
        ! call get_command_argument(1,filename)
        call std_atm_English(0.0,Z_temp,T_temp,P_temp,rho0,a_temp,mu_temp)


        call jsonx_load(filename,j_main)
        ! simulation
        call jsonx_get(j_main, "simulation.time_step[sec]", dt, 0.0)
        call jsonx_get(j_main, "simulation.end_time[sec]", tf, 20.0)
        call jsonx_get(j_main, "simulation.rk4_verbose", rk4_verbose, .false.)
        ! vehicle
        ! thrust
        call jsonx_get(j_main, "vehicle.thrust.Thrust0[lbf]", Thrust0)
        call jsonx_get(j_main, "vehicle.thrust.Ta", Ta)
        ! aerodynamics
        ! reference
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.area[ft^2]", sref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.longitudinal_length[ft]", long_ref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.lateral_length[ft]", lat_ref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.relative_location[ft]", aero_ref_location,0.0,3)
        ! coefficients 
        !CL 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.0", CL0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alpha", CLa)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alphahat", CLahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.qbar", CLqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.elevator", CLde)
        !CS 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.beta", CSb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.pbar", CSpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.alpha_pbar", CSapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rbar", CSrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.aileron", CSda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rudder", CSdr)
        ! CD
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.L0", CDL0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1", CDL1)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1_CL1", CDL2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CS_CS", CDS2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.qbar", CDqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_qbar", CDaqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator", CDde)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_elevator", CDade)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator_elevator", CDde2)
        ! Cl
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.beta", Clb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.pbar", Clpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rbar", Clrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.alpha_rbar", Clarbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.aileron", Clda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rudder", Cldr)
        ! Cm 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.0", Cm0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alpha", Cma)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.qbar", Cmqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alphahat", Cmahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.elevator", Cmde)
        ! Cn
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.beta", Cnb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.pbar", Cnpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_pbar", Cnapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rbar", Cnrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.aileron", Cnda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_aileron", Cnada)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rudder", Cndr)
        ! initial conditions
        call jsonx_get(j_main, "connections.controls.is_use", is_use_controls)
        call mass_inertia() ! computes mass, inertia, and gyroscopic properties of a thing which is weighed at a geopotential altitude of exactly zero. 
        y_init = 0.0
        call jsonx_get(j_main, "initial.time[sec]", t, 0.0)
        call jsonx_get(j_main, "initial.airspeed[ft/sec]", V_initial)
        call jsonx_get(j_main, "initial.altitude[ft]", y_init(9))
        y_init(9) = - y_init(9)
        call jsonx_get(j_main, "initial.Euler_angles[deg]", eul0,0.0,3)
        eul0 = eul0*PI/180.0
        y_init(10:13) = euler_to_quat(eul0)
        call jsonx_get(j_main, "initial.type", init_type)
        is_bank_or_beta_for_shss = "bank"
        if (init_type=="state") then
            call jsonx_get(j_main, "initial.state.angle_of_attack[deg]", alpha_initial)
            alpha_initial = alpha_initial * PI/180.0
            call jsonx_get(j_main, "initial.state.sideslip_angle[deg]", beta_initial)
            beta_initial = beta_initial * PI/180.0
            y_init(1) = V_initial*cos(alpha_initial)*cos(beta_initial)
            y_init(2) = V_initial*sin(beta_initial)
            y_init(3) = V_initial*sin(alpha_initial)*cos(beta_initial)
            call jsonx_get(j_main, "initial.state.angular_rates[deg/s]", angular_rates,0.0,3)
            y_init(4:6) = angular_rates * PI/180.0
            call jsonx_get(j_main, "initial.state.aileron[deg]", controls(1))
            call jsonx_get(j_main, "initial.state.elevator[deg]", controls(2))
            call jsonx_get(j_main, "initial.state.rudder[deg]", controls(3))
            call jsonx_get(j_main, "initial.state.throttle", controls(4))
            controls(1:3) = controls(1:3)*PI/180.0
        else 
            call jsonx_get(j_main, "initial.trim.type", trim_type)
            call jsonx_get(j_main, "initial.trim.elevation_or_climb", is_elevation_or_climb)
            if (is_elevation_or_climb == "climb") then
                call jsonx_get(j_main, "initial.trim.climb_angle[deg]", trim_climb_angle)
                trim_climb_angle = trim_climb_angle*PI/180.0
                trim_elevation_angle = 0.0
            else 
                call jsonx_get(j_main, "initial.trim.elevation_angle[deg]", trim_elevation_angle)
                trim_elevation_angle = trim_elevation_angle*PI/180.0
                trim_climb_angle = 0.0
            end if 
            call jsonx_get(j_main, "initial.trim.solver.finite_difference_step_size", finite_diff_step)
            call jsonx_get(j_main, "initial.trim.solver.relaxation_factor", relax_factor)
            call jsonx_get(j_main, "initial.trim.solver.tolerance", newton_tol)
            call jsonx_get(j_main, "initial.trim.solver.max_iterations", newton_max_iter)
            trim_bank_angle = 0.0
            trim_sideslip_angle = 0.0
            p_wind = 0.0
            trim_azimuth_angle = eul0(3)
            is_trim_sideslip_angle = .false.
            if (trim_type == "sct") then
                call jsonx_get(j_main, "initial.trim.type_sct.bank_angle[deg]", trim_bank_angle)
            else if (trim_type == "shss") then
                call jsonx_get(j_main, "initial.trim.type_shss.bank_or_beta", is_bank_or_beta_for_shss)
                if (is_bank_or_beta_for_shss == "bank") then 
                    call jsonx_get(j_main, "initial.trim.type_shss.bank.bank_angle[deg]", trim_bank_angle)
                else if (is_bank_or_beta_for_shss == "beta") then 
                    is_trim_sideslip_angle = .true.
                    call jsonx_get(j_main, "initial.trim.type_shss.beta.sideslip_angle[deg]", trim_sideslip_angle)
                end if 
            else if (trim_type == "vbr") then 
                call jsonx_get(j_main, "initial.trim.type_vbr.p_wind[deg/s]", p_wind)
                call jsonx_get(j_main, "initial.trim.type_vbr.bank_angle[deg]", trim_bank_angle)
            end if 
            trim_bank_angle = trim_bank_angle*PI/180.0
            trim_sideslip_angle = trim_sideslip_angle*PI/180.0
            p_wind = p_wind*PI/180.0
            ! call jsonx_get(j_main, "initial.trim.elevation_angle[deg]", trim_elevation_angle)
            ! call jsonx_get(j_main, "initial.trim.bank_angle[deg]", trim_bank_angle)
            ! call jsonx_get(j_main, "initial.trim.sideslip_angle[deg]", trim_sideslip_angle)
            alpha_initial = 0.0
            beta_initial = 0.0
            if (trim_type == "shss" .and. is_trim_sideslip_angle) then
                beta_initial = trim_sideslip_angle
            end if 
            write(*,*) 
            write(*,*) 
            write(*,*) "alpha", alpha_initial
            write(*,*) "beta", beta_initial
            write(*,*) "x[ft]", y_init(7)
            write(*,*) "y[ft]", y_init(8)
            write(*,*) "altitude[ft]", y_init(9)
            write(*,*) "trim_type ", trim_type
            write(*,*) "is_bank_or_beta_for_shss ", is_bank_or_beta_for_shss
            write(*,*) "trim_bank_angle", trim_bank_angle
            write(*,*) "trim_elevation_angle", trim_elevation_angle
            write(*,*) "trim_azimuth_angle", trim_azimuth_angle
            write(*,*) 
            write(*,*) 

            write(*,'(a)') 'Trimming Aircraft for '// trim_type
            write(*,'(a,f12.6)') '  --> Azimuth angle set to psi [deg] = ', trim_azimuth_angle*180.0/PI
            write(*,'(a,f12.6)') '  --> Elevation angle set to theta [deg] = ', trim_elevation_angle*180.0/PI
            write(*,'(a,f12.6)') '  --> Bank angle set to phi [deg] = ', trim_bank_angle*180.0/PI

            write(*,'(a,1x,e25.16)') 'Initial theta [deg] = ', trim_elevation_angle*180.0/PI
            write(*,'(a,1x,e25.16)') 'Initial gamma [deg] = ', trim_climb_angle * 180.0/PI
            write(*,'(a,1x,e25.16)') 'Initial phi [deg]   = ', trim_bank_angle*180.0/PI
            write(*,'(a,1x,e25.16)') 'Initial beta [deg]  = ', beta_initial*180.0/PI

            write(*,'(a)') 'Newton Solver Settings:'
            write(*,'(a,1x,e25.16)') 'Finite Difference Step Size = ', finite_diff_step
            write(*,'(a,1x,e25.16)') '          Relaxation Factor = ', relax_factor
            write(*,'(a,1x,e25.16)') '                  Tolerance = ', newton_tol

            trim_array = trim_algorithm(y_init(9), newton_tol)
            if (is_bank_or_beta_for_shss == "beta" .and. trim_type == "shss") then 
                y_init(10:13) = euler_to_quat([trim_array(2), trim_elevation_angle, trim_azimuth_angle])
                write(*,'(A12,1X,E22.13)') "theta[deg]", trim_elevation_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "phi[deg]", trim_array(2)*180.0/PI
                write(*,'(A12,1X,E22.13)') "alpha[deg]", trim_array(1)*180.0/PI
                write(*,'(A12,1X,E22.13)') "beta[deg]", beta_initial*180.0/PI
                write(*,'(A12,1X,E22.13)') "p[deg/s]", trim_array(3)*180.0/PI
                write(*,'(A12,1X,E22.13)') "q[deg/s]", trim_array(4)*180.0/PI
                write(*,'(A12,1X,E22.13)') "r[deg/s]", trim_array(5)*180.0/PI
                write(*,'(A12,1X,E22.13)') "da[deg]", trim_array(6)*180.0/PI
                write(*,'(A12,1X,E22.13)') "de[deg]", trim_array(7)*180.0/PI
                write(*,'(A12,1X,E22.13)') "dr[deg]", trim_array(8)*180.0/PI
                write(*,'(A12,1X,E22.13)') "tau", trim_array(9)
            else
                beta_initial = trim_array(2)
                y_init(10:13) = euler_to_quat([trim_bank_angle, trim_elevation_angle, trim_azimuth_angle])
                write(*,'(A12,1X,E22.13)') "theta[deg]", trim_elevation_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "phi[deg]", trim_bank_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "alpha[deg]", trim_array(1)*180.0/PI
                write(*,'(A12,1X,E22.13)') "beta[deg]", trim_array(2)*180.0/PI
                write(*,'(A12,1X,E22.13)') "p[deg/s]", trim_array(3)*180.0/PI
                write(*,'(A12,1X,E22.13)') "q[deg/s]", trim_array(4)*180.0/PI
                write(*,'(A12,1X,E22.13)') "r[deg/s]", trim_array(5)*180.0/PI
                write(*,'(A12,1X,E22.13)') "da[deg]", trim_array(6)*180.0/PI
                write(*,'(A12,1X,E22.13)') "de[deg]", trim_array(7)*180.0/PI
                write(*,'(A12,1X,E22.13)') "dr[deg]", trim_array(8)*180.0/PI
                write(*,'(A12,1X,E22.13)') "tau", trim_array(9)
                ! write(*,'(A12,1X,E22.13)') "psi[deg]", trim_azimuth_angle*180.0/PI
            end if  
            alpha_initial = trim_array(1)
            y_init(1) = V_initial*cos(alpha_initial)*cos(beta_initial)
            y_init(2) = V_initial*sin(beta_initial)
            y_init(3) = V_initial*sin(alpha_initial)*cos(beta_initial)
            y_init(4:6) = trim_array(3:5)
            controls(1:4) = trim_array(6:9) 
        end if
        ! connections 
        call jsonx_get(j_main, 'connections', j_connections)
        call jsonx_get(j_connections, 'graphics', j_graphics)
        call graphics%init(j_graphics)
        call jsonx_get(j_connections, 'controls', j_controls)
        call connect_controls%init(j_controls)
    end subroutine init

    function trim_algorithm(H_altitude, newton_tol) result(trim_result)
        implicit none
        real, intent(in) :: H_altitude, newton_tol
        real :: trim_result(9)
        real :: alpha, beta, p, q, r, da, de, dr, tau,x,y,z
        real :: pos(3), quat_orientation(4)
        real :: phi, theta, psi
        real :: current_error
        real :: u,v,w, sct_pqr_coeff
        real, allocatable :: DeltaG(:)
        real :: residual(6)
        real :: jacobian(6,6)
        real :: gravity 
        real :: newton_input(6)
        integer :: i, j
        gravity = gravity_English(-H_altitude)
        alpha = 0.0
        allocate(DeltaG(6))
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
            beta = trim_sideslip_angle
        else
            beta = 0.0
        end if 
        u = V_initial*cos(alpha)*cos(beta)
        v = V_initial*sin(beta)
        w = V_initial*sin(alpha)*cos(beta) 
        x = 0.0
        y = 0.0
        z = H_altitude
        da = 0.0
        de = 0.0 
        dr = 0.0
        tau = 0.0
        if (trim_type == "shss") then
            if (is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
            else 
                phi = trim_bank_angle
            end if 
        else
            phi = trim_bank_angle
        end if 
        ! trim_elevation_angle = calc_theta_from_climb_angle(trim_climb_angle, u, v, w, phi)
        if (trim_type == "sct") then 
            sct_pqr_coeff = gravity*sin(trim_bank_angle)*cos(trim_elevation_angle)/&
            (u*cos(trim_elevation_angle)*cos(trim_bank_angle)+w*sin(trim_elevation_angle))
            p = -sct_pqr_coeff*(sin(trim_elevation_angle))
            q = sct_pqr_coeff*(sin(trim_bank_angle)*cos(trim_elevation_angle))
            r = sct_pqr_coeff*(cos(trim_bank_angle)*cos(trim_elevation_angle))
        else if (trim_type == "vbr") then 
            p = (p_wind/V_initial)*u
            q = (p_wind/V_initial)*v
            r = (p_wind/V_initial)*w
        else
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        current_error = 100.0
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
                newton_input = [alpha, phi, da, de, dr, tau]
            else
                newton_input = [alpha, beta, da, de, dr, tau]
        end if
        j = 1
        do while(current_error > newton_tol)
            if (newton_input(6) < 0.0) then 
                newton_input(6) = 0.0
            else if (newton_input(6) > 1.0) then 
                newton_input(6) = 1.0
            end if 
            alpha = newton_input(1)
            if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
                phi = newton_input(2)
            else
                beta = newton_input(2)
            end if 
            ! beta = newton_input(2)
            u = V_initial*cos(alpha)*cos(beta)
            v = V_initial*sin(beta)
            w = V_initial*sin(alpha)*cos(beta) 
            if (is_elevation_or_climb == "climb") then 
                trim_elevation_angle = calc_theta_from_climb_angle(trim_climb_angle, u, v, w, phi)
            end if 
            ! write(*,*) "Trim theta", trim_elevation_angle * 180.0/PI
            if (trim_type == "sct") then  
                sct_pqr_coeff = gravity*sin(trim_bank_angle)*cos(trim_elevation_angle)/&
                (u*cos(trim_elevation_angle)*cos(trim_bank_angle)+w*sin(trim_elevation_angle))
                p = -sct_pqr_coeff*(sin(trim_elevation_angle))
                q = sct_pqr_coeff*(sin(trim_bank_angle)*cos(trim_elevation_angle))
                r = sct_pqr_coeff*(cos(trim_bank_angle)*cos(trim_elevation_angle))
            else if (trim_type == "vbr") then 
                p = (p_wind/V_initial)*u
                q = (p_wind/V_initial)*v
                r = (p_wind/V_initial)*w
            else 
                p = 0.0
                q = 0.0
                r = 0.0
            end if 
            residual = calc_residual(newton_input, p, q, r)
            ! if (rk4_verbose) then
            ! write(*,'(a)') 'Updating rotation rates for ', trim_type
            ! write(*,'(a,1x,e25.16)') 'p [deg/s] = ', (p * 180.0/PI)
            ! write(*,'(a,1x,e25.16)') 'q [deg/s] = ', (q * 180.0/PI)
            ! write(*,'(a,1x,e25.16)') 'r [deg/s] = ', (r * 180.0/PI)
            ! write(*,'(a)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
            ! write(*,'(a,1x,6(e25.16,","))') 'G = ', newton_input
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', residual
            ! write(*,'(a,1x,e25.16)') 'current_error = ', current_error
            ! end if
            ! write(*,*) "residual"
            ! write(*,*) residual
            jacobian = create_jacobian(newton_input, finite_diff_step, p, q, r)
            ! if (rk4_verbose) then
            ! write(*,'(a)') 'Jacobian Matrix ='
            ! do i = 1,6
            !     write(*,'(6(1x,e25.16))') jacobian(i,:)
            ! end do
            ! end if
            ! write(*,*) "jacobian"
            ! write(*,*) jacobian
            call lu_solve(6,jacobian,residual,DeltaG)
            ! write(*,*) "DeltaG"
            ! write(*,*) -DeltaG
            newton_input = newton_input - relax_factor*DeltaG
            ! write(*,*) "New G"
            ! write(*,*) newton_input
            residual = calc_residual(newton_input, p, q, r)
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', residual
            current_error = maxval(abs(residual))
            ! write(*,*) "Iteration, Residual, alpha, beta, p, q, &
            ! r, phi, theta, ail, el, rud, throttle"
            ! write(*,*) j, current_error, newton_input(1)*180.0/PI, &
            ! newton_input(2)*180.0/PI, p*180.0/PI, q*180.0/PI, r*180.0/PI, phi*180.0/PI,trim_elevation_angle*180.0/PI, &
            ! newton_input(3)*180.0/PI, newton_input(4)*180.0/PI, newton_input(5)*180.0/PI, newton_input(6)
            j = j + 1
        end do 
        alpha = newton_input(1)
        beta = newton_input(2)
        da = newton_input(3)
        de = newton_input(4)
        dr = newton_input(5)
        tau = newton_input(6)
        u = V_initial*cos(alpha)*cos(beta)
        v = V_initial*sin(beta)
        w = V_initial*sin(alpha)*cos(beta) 
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
            phi = newton_input(2)
        else 
            phi = trim_bank_angle
        end if
        theta = trim_elevation_angle
        psi = trim_azimuth_angle
        if (trim_type == "sct") then 
            sct_pqr_coeff = gravity*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
            p = -sct_pqr_coeff*(sin(theta))
            q = sct_pqr_coeff*(sin(phi)*cos(theta))
            r = sct_pqr_coeff*(cos(phi)*cos(theta))
        else if (trim_type == "vbr") then 
            p = (p_wind/V_initial)*u
            q = (p_wind/V_initial)*v
            r = (p_wind/V_initial)*w
        else 
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        trim_result = [alpha, beta, p, q, r, da, de, dr, tau]
    end function trim_algorithm

    subroutine mass_inertia() 
        implicit none 
        real :: gravity
        real :: I(3,3)
        real :: det_I
        real :: I_tilde(3,3) 
        real :: a11, a22, a33, a12, a13, a21, a23, a31, a32
       
    end subroutine mass_inertia

    subroutine pseudo_aero(y) 
        implicit none
        real, intent(in) :: y(13)
        real :: da, de, dr, tau
        real :: V, alpha, beta, pbar, qbar, rbar, ahat
        real :: CL1, CL, CS, CD, Cll, Cm, Cn 
        real :: sa, ca, sb, cb
        real :: Z, T, P, rho, a, mu

    end subroutine pseudo_aero

end module sim_m
