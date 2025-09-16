module sim_m
    use adams_m
    implicit none
    ! call test_runge_kutta(t_0, state0)
    contains

    ! subroutine simulation_main(t_0, tf, delta_t, state)
    !     real, intent(in) :: tf, delta_t
    !     real, intent(inout) :: t_0
    !     real, intent(inout), dimension(:) :: state
    !     .while. t_0 < tf do 
    !         state = runge_kutta(t_0, state, delta_t)
    !     end do 
        
    ! end subroutine simulation_main
    subroutine run()
        ! type1, intent(in) :: arg1
        ! type2, intent(out) ::  arg2
    
        
    end subroutine run

    ! test multi-dimensional rk4
    subroutine test_runge_kutta(t_0, state0) 
        implicit none 
        real, intent(in) :: t_0
        real, intent(in), dimension(:) :: state0
        integer :: i
        real :: t, start_time, end_time
        real, dimension(size(state0)) :: state
        call cpu_time(start_time)
        t = t_0
        state = state0
        do i = 1, runge_kutta_array_size
            ! write(*,*) 't, state: ', t, state
            state = runge_kutta(t, state, delta_t)
            t = t + delta_t
        end do
        call cpu_time(end_time)
        write(*,*) 'Rk4 time total [sec]:          ', end_time - start_time
    end subroutine test_runge_kutta

    function differential_equations(t, state) result(res)
        implicit none 
        real, intent(in) :: t
        real, intent(in), dimension(:) :: state
        real :: res(size(state))
        real :: x, z
        x = state(1)
        z = state(2)
        res(1) = t + z/200.0!**2*sin(x)
        res(2) = t + x/200.0!*cos(z)
    end function differential_equations

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
