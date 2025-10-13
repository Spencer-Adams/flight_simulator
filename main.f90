program main
    use sim_m    
    implicit none
    character(100) :: filename 
    call get_command_argument(1, filename)
    call init(filename)
    call run()
end program main