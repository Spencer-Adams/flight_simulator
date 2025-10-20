module sim_m
    use adams_m
    use jsonx_m
    use micro_time_m
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
    real :: V_initial, alpha_initial, beta_initial
    real :: rho0
    logical :: rk4_verbose
    type(json_value), pointer :: j_main
    ! aero coefficients 
    real, allocatable :: aero_ref_location(:)
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
    real, dimension(3) :: eul0 
    integer :: n
    contains
    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        real :: y(13), y1(13)
        real :: Z_temp, P_temp, T_temp, a_temp, mu_temp
        real :: cpu_start_time, cpu_end_time, time1, time2, actual_time, integrated_time
        integer :: io_unit
        logical :: real_time
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        n = 13
        call std_atm_English(0.0,Z_temp,T_temp,P_temp,rho0,a_temp,mu_temp)
        ! initial conditions 
        t = 0.0
        y = y_init

        if (abs(dt)<TOLERANCE) then 
            real_time = .true.
            time1 = get_time()
            y1 = runge_kutta(t, y,dt)
            call quat_norm(y1(10:13))
            time2 = get_time()
            dt = time2-time1 
            write(*,*) "Dt = ", dt
            if (dt == 0.0) then 
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
                write(*,*) "THE DT IS EXACTLY ZERO MY GUY!!!!"
            end if 
            y = y_init 
            t = 0.0
        end if
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
            y1 = runge_kutta(t,y,dt)
            call quat_norm(y1(10:13))
            y = y1
            t = t + dt 
            integrated_time = integrated_time + dt
            write(io_unit,'(14E22.13)') t,y(:)
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

    function cross_product_3D(a, b) result(result)
        implicit none
        real, intent(in) :: a(3), b(3) 
        real :: result(3)
        result(1) = a(2)*b(3) - a(3)*b(2)
        result(2) = a(3)*b(1) - a(1)*b(3)
        result(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product_3D

    function runge_kutta(t_0, state, delta_t) result(state_out)
        implicit none 
        real, intent(in) :: t_0
        real, intent(in), dimension(:) :: state
        real, intent(in) :: delta_t
        real :: state_out(size(state))
        real :: k1(n), k2(n), k3(n), k4(n)
        real :: state_temp(n)
        real :: t_0_plus_delta_t_over_2
        t_0_plus_delta_t_over_2 = t_0 + delta_t_over_2
        k1 = differential_equations(t_0, state)
            state_temp = state
            state_temp = state_temp + delta_t_over_2 * k1
        k2 = differential_equations(t_0_plus_delta_t_over_2, state_temp)
            state_temp = state
            state_temp = state_temp + delta_t_over_2 * k2
        k3 = differential_equations(t_0_plus_delta_t_over_2, state_temp)
            state_temp = state 
            state_temp = state_temp + delta_t * k3 
        k4 = differential_equations(t_0 + delta_t, state_temp)
        state_out = state + delta_t_over_6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    end function runge_kutta

    function differential_equations(time, state) result(res)
        implicit none 
        real, intent(in) :: time
        real, intent(in), dimension(:) :: state
        real :: res(size(state))
        real :: u,v,w,p,q,r,x,y,z
        real :: e0,ex,ey,ez
        real, dimension(3) :: pqr_temp, rot_and_inertia_temp
        real :: gravity
        u = state(1)
        v = state(2)
        w = state(3)
        p = state(4)
        q = state(5)
        r = state(6)
        pqr_temp = [p,q,r]
        x = state(7)
        y = state(8)
        z = state(9)
        e0 = state(10)
        ex = state(11)
        ey = state(12)
        ez = state(13)
        ! write(*,*) "State Vector Coming In"
        ! write(*,*) State
        gravity = gravity_English(-z)
        call pseudo_aero(state)
        rot_and_inertia_temp = & 
        [FM(4) + dot_product(h_gyro(1,:),pqr_temp) + ((Iyyb-Izzb)*q*r-Iyzb*(q**2-r**2)-Ixzb*p*q+Ixyb*p*r)-hdot_gyro(1),&
        FM(5) + dot_product(h_gyro(2,:),pqr_temp) + ((Izzb-Ixxb)*p*r-Ixzb*(r**2-p**2)-Ixyb*q*r+Iyzb*p*q)-hdot_gyro(2),& 
        FM(6) + dot_product(h_gyro(3,:),pqr_temp) + ((Ixxb-Iyyb)*p*q-Ixyb*(p**2-q**2)-Iyzb*p*r+Ixzb*q*r)-hdot_gyro(3)]
        res(1) = 1/mass * FM(1) + gravity * 2 *(ex*ez-ey*e0) + r*v - q*w ! udot body-fixed
        res(2) = 1/mass * FM(2) + gravity * 2 *(ey*ez+ex*e0) + p*w - r*u ! vdot body-fixed
        res(3) = 1/mass * FM(3) + gravity * (ez**2+e0**2-ex**2-ey**2) + q*u - p*v ! wdot body-fixed
        res(4:6) = matmul(Iinv, rot_and_inertia_temp) ! pdot, qdot, rdot body-fixed
        res(7:9) = quat_dependent_to_base((/u,v,w/), (/e0, ex, ey, ez/)) !!! + wind ! xdot ydot zdot earth-fixed
        res(10) = 0.5 * dot_product((/-ex, -ey, -ez/),pqr_temp) !e0
        res(11) = 0.5 * dot_product((/e0, -ez, ey/),pqr_temp) !ex
        res(12) = 0.5 * dot_product((/ez, e0, -ex/),pqr_temp) !ey
        res(13) = 0.5 * dot_product((/-ey, ex, e0/),pqr_temp) !ez
        if (rk4_verbose) then
            write(*,*) "state"
            write(*,*) state
            write(*,*) ""
            write(*,*) "res"
            write(*,'(14E22.13)') res
            write(*,*) ""
        end if
    end function differential_equations

    subroutine init(filename)
        implicit none 
        character(100), intent(in) :: filename
        ! type2, intent(out) ::  arg2
        ! call get_command_argument(1,filename)
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
        y_init = 0.0
        call jsonx_get(j_main, "initial.time[sec]", t, 0.0)
        call jsonx_get(j_main, "initial.airspeed[ft/sec]", V_initial)
        call jsonx_get(j_main, "initial.angle_of_attack[deg]", alpha_initial)
        alpha_initial = alpha_initial * PI/180.0
        call jsonx_get(j_main, "initial.sideslip_angle[deg]", beta_initial)
        beta_initial = beta_initial * PI/180.0
        y_init(1) = V_initial*cos(alpha_initial)*cos(beta_initial)
        y_init(2) = V_initial*sin(beta_initial)
        y_init(3) = V_initial*sin(alpha_initial)*cos(beta_initial)
        call jsonx_get(j_main, "initial.p[deg/s]", y_init(4))
        call jsonx_get(j_main, "initial.q[deg/s]", y_init(5))
        call jsonx_get(j_main, "initial.r[deg/s]", y_init(6))
        y_init(4) = y_init(4) * PI/180.0 
        y_init(5) = y_init(5) * PI/180.0 
        y_init(6) = y_init(6) * PI/180.0
        call jsonx_get(j_main, "initial.altitude[ft]", y_init(9))
        y_init(9) = - y_init(9)

        eul0 = 0.0
        call jsonx_get(j_main, "initial.bank_angle[deg]", eul0(1))
        call jsonx_get(j_main, "initial.elevation_angle[deg]", eul0(2))
        call jsonx_get(j_main, "initial.heading_angle[deg]", eul0(3))
        eul0 = eul0*PI/180.0
        y_init(10:13) = euler_to_quat(eul0)
        call jsonx_get(j_main, "initial.aileron[deg]", controls(1))
        call jsonx_get(j_main, "initial.elevator[deg]", controls(2))
        call jsonx_get(j_main, "initial.rudder[deg]", controls(3))
        call jsonx_get(j_main, "initial.throttle", controls(4))
        controls(1:3) = controls(1:3)*PI/180.0
        call mass_inertia() ! computes mass, inertia, and gyroscopic properties of a thing which is weighed at a geopotential altitude of exactly zero. 
    end subroutine init

    subroutine mass_inertia() 
        implicit none 
        real :: gravity
        real :: I(3,3)
        real :: det_I
        real :: I_tilde(3,3) 
        real :: a11, a22, a33, a12, a13, a21, a23, a31, a32
        I = 0.0
        call jsonx_get(j_main, "vehicle.mass.weight[lbf]", weight)
        call jsonx_get(j_main, "vehicle.mass.Ixx[slug-ft^2]", I(1,1))
        call jsonx_get(j_main, "vehicle.mass.Iyy[slug-ft^2]", I(2,2))
        call jsonx_get(j_main, "vehicle.mass.Izz[slug-ft^2]", I(3,3))
        call jsonx_get(j_main, "vehicle.mass.Ixy[slug-ft^2]", I(1,2))
        call jsonx_get(j_main, "vehicle.mass.Ixz[slug-ft^2]", I(3,1))
        call jsonx_get(j_main, "vehicle.mass.Iyz[slug-ft^2]", I(3,2))
        call jsonx_get(j_main, "vehicle.mass.hx[slug-ft^2/s]", hx) ! in the 3 by 3 h matrix in Eq. 5.4.6
        call jsonx_get(j_main, "vehicle.mass.hy[slug-ft^2/s]", hy)
        call jsonx_get(j_main, "vehicle.mass.hz[slug-ft^2/s]", hz)
        I(1,2) = -I(1,2)
        I(3,1) = -I(3,1)
        I(3,2) = -I(3,2)
        I(2,1) = I(1,2)
        I(1,3) = I(3,1)
        I(2,3) = I(3,2)
        Ixxb = I(1,1)
        Iyyb = I(2,2)
        Izzb = I(3,3)
        Ixyb = I(1,2)
        Ixzb = I(1,3)
        Iyzb = I(2,3)
        a11 = I(1,1)
        a22 = I(2,2)
        a33 = I(3,3)
        a12 = I(1,2)
        a13 = I(1,3)
        a21 = I(2,1)
        a23 = I(2,3)
        a31 = I(3,1)
        a32 = I(3,2)
        gravity = gravity_English(0.0)
        mass = weight/gravity
        det_I = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 &
        - a13*a22*a31 - a12*a21*a33 - a11*a23*a32
        I_tilde(1,1) = a22*a33 - a23*a32
        I_tilde(1,2) = a13*a32 - a12*a33
        I_tilde(1,3) = a12*a23 - a13*a22
        I_tilde(2,1) = a23*a31 - a21*a33
        I_tilde(2,2) = a11*a33 - a13*a31
        I_tilde(2,3) = a13*a21 - a11*a23
        I_tilde(3,1) = a21*a32 - a22*a31
        I_tilde(3,2) = a12*a31 - a11*a32
        I_tilde(3,3) = a11*a22 - a12*a21
        Iinv = 1/(det_I)*I_tilde
        h_gyro = 0.0
        h_gyro(1,2) = -hz
        h_gyro(1,3) = hy
        h_gyro(2,1) = hz
        h_gyro(2,3) = -hx
        h_gyro(3,1) = -hy
        h_gyro(3,2) = hx
        hdot_gyro = 0.0
    end subroutine mass_inertia

    subroutine pseudo_aero(y) 
        implicit none
        real, intent(in) :: y(13)
        real :: da, de, dr, tau
        real :: V, alpha, beta, pbar, qbar, rbar, ahat
        real :: CL1, CL, CS, CD, Cll, Cm, Cn 
        real :: sa, ca, sb, cb
        real :: Z, T, P, rho, a, mu
        ahat = 0.0

        da = controls(1)
        de = controls(2)
        dr = controls(3)
        tau = controls(4)

        call std_atm_English(-y(9), Z, T, P, rho, a, mu)

        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        alpha = atan2(y(3), y(1)) ! Eq. 3.4.4
        beta = asin(y(2)/V) ! Eq. 3.4.5
        pbar = 0.5*y(4)*lat_ref/(V)
        qbar = 0.5*y(5)*long_ref/(V)
        rbar = 0.5*y(6)*lat_ref/(V)

        CL1 = CL0 +CLa*alpha
        CL = CL1 + CLqbar*qbar+CLahat*ahat + CLde*de
        CS = CSb*beta + (CSpbar+CSapbar*alpha)*pbar + CSrbar*rbar + CSda*da + CSdr*dr
        CD = CDL0 + CDL1*CL1 + CDL2*CL1**2 + CDS2*CS**2 + (CDqbar + CDaqbar*alpha)*qbar + (CDde + CDade*alpha)*de + CDde2*de**2
        Cll = Clb*beta + Clpbar*pbar + (Clrbar + Clarbar*alpha)*rbar + Clda*da + Cldr*dr
        Cm = Cm0 + Cma*alpha + Cmqbar*qbar + Cmahat*ahat + Cmde*de 
        Cn = Cnb*beta + (Cnpbar + Cnapbar*alpha)*pbar + Cnrbar*rbar + (Cnda + Cnada*alpha)*da + Cndr*dr
        sa = sin(alpha)
        ca = cos(alpha)
        sb = sin(beta)
        cb = cos(beta)
        FM(1) = CL*sa-CS*ca*sb-CD*ca*cb
        FM(2) = CS*cb-CD*sb
        FM(3) = -CL*ca-CS*sa*sb-CD*sa*cb
        FM(4) = lat_ref*Cll
        FM(5) = long_ref*Cm
        FM(6) = lat_ref*Cn
        FM = 0.5*rho*V**2*sref*FM
        FM(1) = FM(1) + tau*Thrust0*(rho/rho0)**Ta
        FM(4:6) = FM(4:6) + cross_product_3D(aero_ref_location, FM(1:3)) ! body-fixed axis
        if (rk4_verbose) then
            write(*,*) "FM"
            write(*,*) FM
        end if
    end subroutine pseudo_aero

end module sim_m
