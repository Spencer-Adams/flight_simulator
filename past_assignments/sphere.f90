module sim_m
    use adams_m
    implicit none
    ! Variables within sim_m
    real :: mass
    real :: I(3,3)
    real :: Ixxb, Iyyb, Izzb, Ixyb, Ixzb, Iyzb
    real :: Iinv(3,3)
    real :: h_gyro(3,3)
    real :: hdot_gyro(3)
    real :: FM(6)
    real :: weight 
    real :: r2 
    real :: thick

    contains

    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        real :: t, dt, tf, y(13), y1(13),eul(3)
        integer :: io_unit
        ! thick = 0.00131 ! ft
        ! r2 = 0.13084/2.0 ! ft
        ! weight = 0.006 ! lbf 
        thick = 0.008333333333333333 ! ft
        r2 = 0.360833333333333 ! ft
        weight = 0.99208018236122 ! lbf

        call mass_inertia(mass, I, Iinv, h_gyro, hdot_gyro, 0.0, weight, r2, thick) ! placeholder 
        Ixxb = I(1,1)
        Iyyb = I(2,2)
        Izzb = I(3,3)
        Ixyb = -I(1,2)
        Ixzb = -I(1,3)
        Iyzb = -I(2,3)

        t = t_0
        dt = delta_t
        tf = t_stop

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
        y(9) = z_0 ! ft alt
        eul = eul0
        y(10:13) = euler_to_quat(eul)
        open(newunit=io_unit, file='sphere_output.txt', status='replace', action='write')
        call quat_norm(y(10:13))
        write(io_unit,*) "        t[s]                   u[ft/s]               v[ft/s]"// &
                   "                 w[ft/s]                p[rad/s]             q[rad/s]"  // &
                   "            r[rad/s]                x[ft]                  y[ft]" // &                   
                   "                  z[ft]                    e0                   ex" // &
                   "                      ey                   ez"
        write(io_unit,'(14E22.13)') t,y(:) ! 14E17.9 is for 14 numbers in scientific notation, each 17 characters wide with 9 after the decimal
        do while(y(9)<0.0 .or. t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
            y1 = runge_kutta(t,y,dt)
            call quat_norm(y1(10:13))
            y = y1
            t = t + dt 
            write(io_unit,'(14E22.13)') t,y(:)
        end do
    end subroutine run

    subroutine mass_inertia(mass, I, Iinv, h_gyro, hdot_gyro, H, weight, r2,thick)
        implicit none 
        real, intent(in) :: weight, r2, thick, H
        real, intent(inout) ::  mass
        real, intent(inout) :: I(3,3), Iinv(3,3), h_gyro(3,3) 
        real, intent(inout) :: hdot_gyro(3)
        real :: gravity, r1, coeff

        gravity = gravity_English(H)
        mass = weight/gravity
        r1 = r2-thick
        coeff = mass * (2*(r2**5-r1**5))/(5*(r2**3-r1**3))
        I = coeff*identity
        Iinv = (1.0/coeff) *identity
        h_gyro = 0.0
        hdot_gyro = 0.0
    end subroutine mass_inertia

    function differential_equations(time, state) result(res)
        implicit none 
        real, intent(in) :: time
        real, intent(in), dimension(:) :: state
        real :: res(size(state))
        real :: u_,v_,w_,p_,q_,r_,x_,y_,z_
        real :: e0_,ex_,ey_,ez_
        real, dimension(3) :: M_c, F_c, V_c, pqr_temp, rot_and_inertia_temp, quat_body_to_earth_temp
        real :: Z, T, P, rho, a, mu, gravity, t00
        t00 = time 
        u_ = state(1)
        v_ = state(2)
        w_ = state(3)
        p_ = state(4)
        q_ = state(5)
        r_ = state(6)
        pqr_temp = (/p_,q_,r_/)
        x_ = state(7)
        y_ = state(8)
        z_ = state(9)
        e0_ = state(10)
        ex_ = state(11)
        ey_ = state(12)
        ez_ = state(13)
        quat_body_to_earth_temp = quat_dependent_to_base((/u_,v_,w_/), (/e0_, ex_, ey_, ez_/))
        V_c(1) = u_
        V_c(2) = v_ 
        V_c(3) = w_
        gravity = gravity_English(-z_)
        call std_atm_English(-z_, Z, T, P, rho, a, mu)
        call pseudo_aero(M_c, F_c, rho, mu, r2, V_c)
        rot_and_inertia_temp = & 
        (/M_c(1) + dot_product(h_gyro(1,:),pqr_temp) + ((Iyyb-Izzb)*q_*r_+Iyzb*(q_**2-r_**2)+Ixzb*p_*q_-Ixyb*p_*r_)-hdot_gyro(1),&
        M_c(2) + dot_product(h_gyro(2,:),pqr_temp) + ((Izzb-Ixxb)*p_*r_+Ixzb*(r_**2-p_**2)+Ixyb*q_*r_-Iyzb*p_*q_)-hdot_gyro(2),& 
        M_c(3) + dot_product(h_gyro(3,:),pqr_temp) + ((Ixxb-Iyyb)*p_*q_+Ixyb*(p_**2-q_**2)+Iyzb*p_*r_-Ixzb*q_*r_)-hdot_gyro(3)/)
        res(1) = 1/mass * F_c(1) + gravity * 2 *(ex_*ez_-ey_*e0_) + r_*v_ - q_*w_
        res(2) = 1/mass * F_c(2) + gravity * 2 *(ey_*ez_+ex_*e0_) + p_*w_ - r_*u_
        res(3) = 1/mass * F_c(3) + gravity * (ez_**2+e0_**2-ex_**2-ey_**2) + q_*u_ - p_*v_
        res(4:6) = matmul(Iinv, rot_and_inertia_temp)
        res(7) = quat_body_to_earth_temp(1) !!! + wind
        res(8) = quat_body_to_earth_temp(2) !!! + wind
        res(9) = quat_body_to_earth_temp(3) !!! + wind
        res(10) = 0.5 * dot_product((/-ex_, -ey_, -ez_/),pqr_temp)
        res(11) = 0.5 * dot_product((/e0_, -ez_, ey_/),pqr_temp)
        res(12) = 0.5 * dot_product((/ez_, e0_, -ex_/),pqr_temp)
        res(13) = 0.5 * dot_product((/-ey_, ex_, e0_/),pqr_temp)
    end function differential_equations

    subroutine pseudo_aero(M_c, F_c, density, mu, r2, V_c)
        implicit none
        real, intent(in) :: density, mu, r2
        real, intent(in), dimension(3) :: V_c
        real, intent(inout), dimension(3) ::  M_c, F_c
        real :: V_mag,  C_D, Reyn
        real, dimension (3) :: unit_V
        M_c = 0.0
        V_mag = sqrt(V_c(1)**2 + V_c(2)**2 + V_c(3)**2)
        unit_V = V_c/V_mag
        Reyn = 2*density*V_mag*r2/mu
        if (Reyn < 0.01) then 
            C_D = 2405.0
        else if (0.01 <= Reyn .and. Reyn <= 450000.0) then 
            C_D = 24/Reyn + 6/(1+sqrt(Reyn)) + 0.4
        else if (450000.0 <= Reyn .and. Reyn <= 560000.0) then 
            C_D = 1.0*10.0**29*Reyn**(-5.211)
        else if (560000.0 <= Reyn .and. Reyn <= 14000000.0) then 
            C_D = -2.0*10.0**(-23)*Reyn**3 - 1.0*10.0**(-16)*Reyn**2 + 9.0*10.0**(-9)*Reyn + 0.069
        else 
            C_D = 0.12
        end if 
        F_c = -0.5*density*V_mag**2*PI*r2**2*C_D*unit_V
    end subroutine pseudo_aero

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

end module sim_m
