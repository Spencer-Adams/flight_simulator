module sim_m
    use adams_m
    use jsonx_m
    implicit none
    ! Variables within sim_m
    real :: mass
    real :: I(3,3)
    real :: Ixxb, Iyyb, Izzb, Ixyb, Ixzb, Iyzb
    real :: Iinv(3,3)
    real :: h_gyro(3,3)
    real :: hdot_gyro(3)
    real :: FM(6)
    real :: weight, long_ref_length, lat_ref_length, ref_area !! might switch this out for more general stuff.
    real ::  t, u_0, v_0, w_0, p_0, q_0, r_0, x_0, y_0, z_0
    real :: bank_angle, elevation_angle, azimuth_angle
    real :: dt, tf, delta_t_over_2, delta_t_over_6
    real :: CL_alpha, CD_0, CD_2, C_malpha, C_mqbar, C_l0, C_lpbar
    real, dimension(3) :: eul0 
    integer :: n
    contains

    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        type(json_value), pointer :: j_main
        real :: y(13), y1(13),eul(3)
        integer :: io_unit
        character(100) :: filename
        logical :: is_straight_fletchings 
        !!! stuff pertaining to arrow only !!!
        CL_alpha = 4.929
        CD_0 = 5.096
        CD_2 = 48.138
        C_malpha = -2.605
        C_mqbar = -9.06
        C_l0 = 3.223
        C_lpbar = -5.378
        call get_command_argument(1,filename)
        call jsonx_load(filename,j_main)
        call jsonx_get(j_main, "simulation.time_step[sec]", dt, 0.01)
        call jsonx_get(j_main, "simulation.weight[lbf]", weight, 0.01)
        call jsonx_get(j_main, "simulation.long_ref_length[ft]", long_ref_length, 2.3)
        call jsonx_get(j_main, "simulation.lat_ref_length[ft]", lat_ref_length, 2.3)
        call jsonx_get(j_main, "simulation.ref_area[ft2]", ref_area, 0.000218)
        call jsonx_get(j_main, "simulation.final_time[sec]", tf, 20.0)
        call jsonx_get(j_main, "initial.is_straight_fletchings", is_straight_fletchings, .true.)
        call jsonx_get(j_main, "initial.time[sec]", t, 0.0)
        call jsonx_get(j_main, "initial.u_body_fixed[ft/s]", u_0, 210.0)
        call jsonx_get(j_main, "initial.v_body_fixed[ft/s]", v_0, 0.0)
        call jsonx_get(j_main, "initial.w_body_fixed[ft/s]", w_0, 0.0)
        call jsonx_get(j_main, "initial.p_body_fixed[deg/s]", p_0, 0.0)
        call jsonx_get(j_main, "initial.q_body_fixed[deg/s]", q_0, 0.0)
        call jsonx_get(j_main, "initial.r_body_fixed[deg/s]", r_0, 0.0)
        call jsonx_get(j_main, "initial.x_earth_fixed[ft]", x_0, 0.0)
        call jsonx_get(j_main, "initial.y_earth_fixed[ft]", y_0, 0.0)
        call jsonx_get(j_main, "initial.z_earth_fixed[ft]", z_0, 20.0)
        call jsonx_get(j_main, "initial.bank_angle[deg]", bank_angle, 0.0)
        call jsonx_get(j_main, "initial.elevation_angle[deg]", elevation_angle, 0.0)
        call jsonx_get(j_main, "initial.azimuth_angle[deg]", azimuth_angle, 0.0)
        eul0 = [bank_angle, elevation_angle, azimuth_angle]
        eul0(1) = eul0(1) * PI/180.0 
        eul0(2) = eul0(2) * PI/180.0 
        eul0(3) = eul0(3) * PI/180.0 
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        n = 13

        CL_alpha = 4.929
        CD_0 = 5.096
        CD_2 = 48.138
        C_malpha = -2.605
        C_mqbar = -9.06
        if (is_straight_fletchings) then 
            C_l0 = 0.0
        else
            C_l0 = 3.223
        end if
        C_lpbar = -5.378

        call mass_inertia(mass, I, Iinv, h_gyro, hdot_gyro, 0.0, weight) ! placeholder 
        Ixxb = I(1,1)
        Iyyb = I(2,2)
        Izzb = I(3,3)
        Ixyb = -I(1,2)
        Ixzb = -I(1,3)
        Iyzb = -I(2,3)

        ! initial conditions 
        y = 0.0
        y(1) = u_0 ! ft/s
        y(2) = v_0 ! ft/s
        y(3) = w_0 ! ft/s
        y(4) = p_0 ! rad/s
        y(5) = q_0 ! rad/s
        y(6) = r_0 ! rad/s
        y(7) = x_0 ! ft
        y(8) = y_0 ! ft
        y(9) = -z_0 ! ft alt
        eul = eul0
        y(10:13) = euler_to_quat(eul)
        open(newunit=io_unit, file='arrow_output.txt', status='replace', action='write')
        call quat_norm(y(10:13))
        write(io_unit,*) "        t[s]                   u[ft/s]               v[ft/s]"// &
                   "                 w[ft/s]                p[rad/s]             q[rad/s]"  // &
                   "            r[rad/s]                x[ft]                  y[ft]" // &                   
                   "                  z[ft]                    e0                   ex" // &
                   "                      ey                   ez"
        write(io_unit,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        ! write(*,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        do while(y(9)<0.0 .and. t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
            y1 = runge_kutta(t,y,dt)
            call quat_norm(y1(10:13))
            y = y1
            t = t + dt 
            write(io_unit,'(14E22.13)') t,y(:)
            ! write(*,'(14E22.13)') t,y(:)
        end do
    end subroutine run

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
        real :: u_,v_,w_,p_,q_,r_,x_,y_,z_
        real :: e0_,ex_,ey_,ez_
        real, dimension(3) :: M_c, F_c, V_c, pqr_temp, rot_and_inertia_temp
        real :: Z, T, P, rho, a, mu, gravity, t00
        t00 = time 
        ! write(*,*) "state"
        ! write(*,*) state
        ! write(*,*) ""
        u_ = state(1)
        v_ = state(2)
        w_ = state(3)
        p_ = state(4)
        q_ = state(5)
        r_ = state(6)
        pqr_temp = (/ p_,q_,r_ /)
        x_ = state(7)
        y_ = state(8)
        z_ = state(9)
        e0_ = state(10)
        ex_ = state(11)
        ey_ = state(12)
        ez_ = state(13)
        V_c(1) = u_
        V_c(2) = v_ 
        V_c(3) = w_
        ! write(*,*) "State Vector Coming In"
        ! write(*,*) State
        gravity = gravity_English(-z_)
        call std_atm_English(-z_, Z, T, P, rho, a, mu)
        call pseudo_aero(M_c, F_c, rho, V_c, p_, q_, r_)
        ! write(*,*) "fc, mc"
        ! write(*,*) f_c, m_c
        ! write(*,*) ""
        rot_and_inertia_temp = & 
        (/ M_c(1) + dot_product(h_gyro(1,:),pqr_temp) + ((Iyyb-Izzb)*q_*r_+Iyzb*(q_**2-r_**2)+Ixzb*p_*q_-Ixyb*p_*r_)-hdot_gyro(1),&
        M_c(2) + dot_product(h_gyro(2,:),pqr_temp) + ((Izzb-Ixxb)*p_*r_+Ixzb*(r_**2-p_**2)+Ixyb*q_*r_-Iyzb*p_*q_)-hdot_gyro(2),& 
        M_c(3) + dot_product(h_gyro(3,:),pqr_temp) + ((Ixxb-Iyyb)*p_*q_+Ixyb*(p_**2-q_**2)+Iyzb*p_*r_-Ixzb*q_*r_)-hdot_gyro(3) /)
        res(1) = 1/mass * F_c(1) + gravity * 2 *(ex_*ez_-ey_*e0_) + r_*v_ - q_*w_ ! udot body-fixed
        res(2) = 1/mass * F_c(2) + gravity * 2 *(ey_*ez_+ex_*e0_) + p_*w_ - r_*u_ ! vdot body-fixed
        res(3) = 1/mass * F_c(3) + gravity * (ez_**2+e0_**2-ex_**2-ey_**2) + q_*u_ - p_*v_ ! wdot body-fixed
        res(4:6) = matmul(Iinv, rot_and_inertia_temp) ! pdot, qdot, rdot body-fixed
        res(7:9) = quat_dependent_to_base((/u_,v_,w_/), (/e0_, ex_, ey_, ez_/)) !!! + wind ! xdot ydot zdot earth-fixed
        res(10) = 0.5 * dot_product((/-ex_, -ey_, -ez_/),pqr_temp) !e0
        res(11) = 0.5 * dot_product((/e0_, -ez_, ey_/),pqr_temp) !ex
        res(12) = 0.5 * dot_product((/ez_, e0_, -ex_/),pqr_temp) !ey
        res(13) = 0.5 * dot_product((/-ey_, ex_, e0_/),pqr_temp) !ez
        ! write(*,*) "res"
        ! write(*,'(14E22.13)') res
        ! write(*,*) ""
    end function differential_equations

    subroutine mass_inertia(mass, I, Iinv, h_gyro, hdot_gyro, H, weight) !!!! need to change mass and inertia to be arrow based
        implicit none 
        real, intent(in) :: weight, H
        real, intent(inout) ::  mass
        real, intent(inout) :: I(3,3), Iinv(3,3), h_gyro(3,3) 
        real, intent(inout) :: hdot_gyro(3)
        real :: gravity
        gravity = gravity_English(H)
        mass = weight/gravity
        I = 0.0
        Iinv = 0.0
        I(1,1) = 0.0000194 ! slug*ft^2
        I(2,2) = 0.00097 ! slug*ft^2
        I(3,3) = 0.00097 ! slug*ft^2
        Iinv(1,1) = 1/I(1,1)
        Iinv(2,2) = 1/I(2,2)
        Iinv(3,3) = 1/I(3,3)
        h_gyro = 0.0
        hdot_gyro = 0.0
    end subroutine mass_inertia

    subroutine pseudo_aero(M_c, F_c, density, V_c, p, q, r) !!!! need to change to arrow based
        implicit none
        real, intent(in) :: density, p, q, r
        real, intent(in), dimension(3) :: V_c
        real, intent(inout), dimension(3) ::  M_c, F_c
        real :: V_mag, alpha, beta, beta_f
        real :: pbar, qbar, rbar, CL, CS, CD, C_l, C_m, C_n
        real :: dynamic ! 0.5*rho*V^2*Sw
        alpha = atan2(V_c(3), V_c(1)) ! Eq. 3.4.4
        V_mag = sqrt(V_c(1)**2 + V_c(2)**2 + V_c(3)**2)
        beta = asin(V_c(2)/V_mag) ! Eq. 3.4.5
        beta_f = atan2(V_c(2), V_c(1)) ! Eq. 3.4.13
        pbar = p*lat_ref_length/(2*V_mag)
        qbar = q*long_ref_length/(2*V_mag)
        rbar = r*lat_ref_length/(2*V_mag)
        CL = CL_alpha*alpha 
        CS = -CL_alpha*beta_f 
        CD = CD_0 + CD_2*CL**2 + CD_2*CS**2
        C_l = C_l0 + C_lpbar*pbar
        C_m = C_malpha*alpha + C_mqbar*qbar 
        C_n = -C_malpha*beta_f + C_mqbar*rbar 
        dynamic =  0.5*density*V_mag**2*ref_area
        ! Eq. 5.2.3 in book. 
        F_c(1) = dynamic*(CL*sin(alpha)-CS*cos(alpha)*sin(beta)-CD*cos(alpha)*cos(beta))
        F_c(2) = dynamic*(CS*cos(beta)-CD*sin(beta))
        F_c(3) = dynamic*(-CL*cos(alpha)-CS*sin(alpha)*sin(beta)-CD*sin(alpha)*cos(beta))
        ! Eq. 5.2.4 in book
        M_c(1) = dynamic*lat_ref_length*C_l
        M_c(2) = dynamic*long_ref_length*C_m
        M_c(3) = dynamic*lat_ref_length*C_n

    end subroutine pseudo_aero

end module sim_m
