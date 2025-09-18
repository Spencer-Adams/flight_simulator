program main
    use adams_m
    implicit none
    integer :: i,j
    real :: Z,T,P,rho,a,mu,g
    real :: start_time, end_time
    integer :: io_unit
    ! real :: R_EZ_eng, R_E_eng
    ! R_EZ_eng = m_to_ft(R_EZ)
    ! write(*,*) "English R_EZ: ", R_EZ_eng

    ! R_E_eng = ft_to_m(R_E_english)
    ! write(*,*) " R_E: ", R_E_eng


    open(newunit=io_unit, file='3.13.7_output.txt', status='replace', action='write')
    write(io_unit,*) 'SI atmospheric test'
    call cpu_time(start_time)
    do i=0,90000,100
        do j = 0,10000
            call std_atm_SI(real(i),Z,T,P,rho,a,mu)
            g = gravity_SI(real(i))
        end do
        write(io_unit,'(7ES25.11)') real(i),Z,T,P,rho,a,mu,g ! the 7ES25.11 is for 7 scientific notation numbers, each 25 characters wide with 11 after the decimal
    end do
    call cpu_time(end_time)
    print *, 'SI time total [sec]:          ', end_time - start_time

    write(io_unit,*) 'English atmospheric test'
    call cpu_time(start_time)
    do i=0,200000,200
        do j=1,10000
            call std_atm_English(real(i),Z,T,P,rho,a,mu)
            g = gravity_English(real(i))
        end do
        write(io_unit,'(7ES25.11)') real(i),Z,T,P,rho,a,mu,g
    end do
    call cpu_time(end_time)
    print *, 'English time total [sec]:     ', end_time - start_time

    close(io_unit)
end program main
