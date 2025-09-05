program function_fun 
    implicit none 

    real, dimension(:), allocatable :: my_array
    integer :: N, i

    ! initialize 
    N = 100
    ! allocate(my_array(100))
    my_array = generate_array(100)
    ! process 
    call process_array(my_array)
    call process_array(my_array)
    call process_array(my_array)
    
    write(*,*) my_array

    contains

    subroutine process_array(array)
        ! processes the array
        implicit none 
        real,dimension(:),allocatable,intent(inout) :: array 
        integer :: j
        !process the array 
        do j = 1, size(array)
            array(j) = 2.*array(j) + 5
        end do 
    end subroutine process_array

    function generate_array(num_elements) result(new_array)
        implicit none 

        integer, intent(in) :: num_elements
        real,dimension(:),allocatable :: new_array
        ! Generate
        !process the array 
        allocate(new_array(num_elements))
        do i = 1, num_elements
            new_array(i) = i
        end do 

    end function generate_array

end program function_fun

! program function_fun 
!     implicit none 

!     real, dimension(:), allocatable :: my_array
!     integer :: N, i

!     ! initialize 
!     N = 100

!     allocate(my_array(100))
!     do i =1, N
!       my_array(i) = i
!     end do 
!     ! process 
!     do i =1, N
!         my_array(i) = 2.*my_array(i) + 5.
!     end do 
!     ! process again
!     do i = 1,N 
!         my_array(i) = 2.*my_array(i) + 5.
!     end do 
    
! end program function_fun

! program helper
!     implicit none
!     real, dimension(:,:), allocatable :: my_array, my_other_array
!     integer :: i   ! <--- declare loop/index variable
!     ! Allocate with 5 rows and 8 columns
!     allocate(my_array(4,4))
!     ! Fill with numbers 12, 13, 14, ...
!     my_array = reshape([(i, i=1, 5*8-1)], shape(my_array))
!     ! Allocate copy
!     allocate(my_other_array, source=my_array)

!     ! Print results
!     write(*,*) "my_array:"
!     write(*,'(8F6.1)') my_array   ! print 8 columns wide

!     write(*,*) "Size of my_array:", size(my_array)

!     write(*,*) "my_other_array:"
!     write(*,'(8F6.1)') my_other_array

!     write(*,*) "Size of my_other_array:", size(my_other_array)
! end program helper

! ! allocatable arrays
! program helper
!     ! Prints hello_world
!     implicit none
!     real, dimension(:,:), allocatable :: my_array, my_other_array
!     ! Initialize variables
!     my_array = 1
!     allocate(my_array(2,2), source=12.)
!     allocate(my_other_array, source = my_array) ! if both arrays are the same size (n,m) then this works, if not, this will NOT work
!     write(*,*) 
!     write(*,*) my_array
!     write(*,*) size(my_array)
!     write(*,*) 
!     write(*,*) my_other_array
!     write(*,*) size(my_other_array)
! end program helper

! basic array - dimension(5,5) means five rows, five columns  
! program helper
!     ! Prints hello_world
!     implicit none
!     integer, dimension(5,5) :: my_array
!     ! Initialize variables
!     my_array = 1
!     write(*,*) 
!     write(*,*) my_array
!     write(*,*) size(my_array)
! end program helper

!!!!!!!DIFFERENT TYPES!!!!!!!!!
! program helper
!     ! Prints hello_world
!     implicit none
!     integer :: my_int
!     real :: my_real 
!     character :: my_char
!     logical :: my_logical
!     complex :: my_complex 
!     ! Initialize variables
!     my_int = 1
!     my_real = 10.
!     my_char = 'a'
!     my_logical = .true.
!     my_complex = (1., 1.)
!     write(*,*) 
!     write(*,*) my_int
!     write(*,*) my_real
!     write(*,*) my_char
!     write(*,*) my_logical
!     write(*,*) my_complex
!     write(*,*) my_complex*my_real
! end program helper

! program helper
!     ! Prints hello_world
!     implicit none
!     real :: x
!     x = 1.
!     write(*,*) "Hello World!", x

! end program helper
