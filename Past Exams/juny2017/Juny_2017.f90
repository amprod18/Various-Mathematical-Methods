program juny2017
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: i, j, k, N(2), infected, n_points
    integer, allocatable :: people(:)
    real(dp) :: t_max, p(3), step, dprevalence
    real(dp), allocatable :: t_RK4(:), prev_RK4(:)
    external :: dprevalence

    t_max = 40._dp
    N = [10000, 100000]
    p = [2._dp/3._dp, 1._dp/2._dp, 2._dp/5._dp]

    n_points = 400
    step = t_max/real(n_points - 1, dp)
    allocate(t_RK4(n_points), prev_RK4(n_points))

    open(1, file='res.dat')
    
    do i = 1, size(N)
        allocate(people(N(i)))
        do j = 1, size(p)
            people = 1
            infected = N(i)
            call simulate_epidemic(N(i), infected, t_max, people, p(j))


            write(1, *) '# RK4 with N=', N, 'and p=', p
            do k = 1, n_points
                t_RK4(k) = step*(k-1)
            end do
            call RungeKutta4order(n_points, step, 1._dp, prev_RK4, p(j), dprevalence)
            do k = 1, n_points
                write(1, '(2f20.14)') t_RK4(k), prev_RK4(k)
            end do
            write(1, *) ' '
            write(1, *) ' '
        end do
        deallocate(people)
    end do

    close(1)

    call system('gnuplot -p figE1.gnu')
    call system('gnuplot -p figE2.gnu')
    call system('gnuplot -p figE3.gnu')
end program juny2017

subroutine simulate_epidemic(N, infected, t_max, people, p)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: N, i, index, index2, status, infected, people(N)
    real(dp) :: x, t_max, p, step, t, step_c, prev
    write(1, *) '# Simulation with N=', N, 'and p=', p
    t = 0._dp
    status = 0
    step_c = 0._dp
    i = 0
    do while (t .le. t_max)
        i = i+1
        ! print *, people
        do while (status .ne. 1)
            call random_number(x)
            index = floor(x*(N-1) + 1)
            status = people(index)
        end do
        index2 = index
        status = 0
        call random_number(x)
        if (x .le. p) then
            people(index) = 0
            infected = infected - 1
        else
            do while (index2 .eq. index)
                call random_number(x)
                index2 = floor(x*(N-1) + 1)
            end do
            if (people(index2) .eq. 0) then
                people(index2) = 1
                infected = infected + 1
            end if
        end if 
        if (infected .lt. 2) then
            write(1, '(2f20.14)') t, 0._dp
            exit
        end if
        step = p/real(infected, dp)
        step_c = step_c + step
        t = t + step
        prev = real(infected, dp)/real(N, dp)
        
        if (step_c .ge. 0.01_dp) then
            write(1, '(2f20.14)') t, prev
            step_c = 0._dp
        end if
        
    end do
    write(1, *) ' '
    write(1, *) ' '

    print *, 'Ended Pandemic in', i, ' iterations'
end subroutine simulate_epidemic

subroutine update_ki(nequs, x, y, ki, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: nequs
    real(dp) :: x, y, func, ki(nequs)

    ki(nequs) = func(x, y)
end subroutine update_ki

subroutine RungeKutta4order(n_points, step, y0, y, p, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i
    real(dp):: step, k1, k2, k3, k4, y0, y(n_points), func, p
    external :: func
    
    ! --- RK4 ---
    ! k1 = f(x0, y0)
    ! k2 = f(x0 + h/2, y0 + h·k1/2)
    ! k3 = f(x0 + h/2, y0 + h·k2/2)
    ! k4 = f(x0 + h, y0 + h·k3)
    ! y1 = y0 + h(k1 + 2·k2 + 2·k3 + k4)/6
    ! ------------

    k1 = func(y0, p)
    k2 = func(y0*k1/2._dp, p)
    k3 = func(y0*k2/2._dp, p)
    k4 = func(y0*k3, p)
    y(1) = y0 + step*(k1 + 2._dp*k2 + 2._dp*k3 + k4)/6._dp

    do i = 2, n_points
        k1 = func(y(i - 1), p)
        k2 = func(y(i - 1)*k1/2._dp, p)
        k3 = func(y(i - 1)*k2/2._dp, p)
        k4 = func(y(i - 1)*k3, p)
        y(i) = y(i - 1) + step*(k1 + 2._dp*k2 + 2._dp*k3 + k4)/6._dp
    end do
end subroutine RungeKutta4order

function dprevalence(prev, p) result(dprev)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp):: prev, dprev, p

    dprev = -prev + (1._dp-p)*prev*(1._dp-prev)/p
end function dprevalence