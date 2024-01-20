program Derivatives
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: xmin ,xmax, func, dfunc, ddfunc
    real(dp), allocatable :: x(:), fx(:), dfx(:), ddfx(:), df_for(:), df_back(:), df_cen(:), df_5p(:), ddf_num(:)

    ! ------------------------ Data Generation ------------------------
    n_points = 20
    xmin = 0._dp
    xmax = 3._dp
    allocate(x(n_points), fx(n_points), dfx(n_points), ddfx(n_points), df_for(n_points), df_back(n_points), &
             df_cen(n_points), df_5p(n_points), ddf_num(n_points))

    open(1, file='Derivatives_res.dat')

    write(1, *) '# Analytic Function'
    do i = 1, n_points
        x(i) = xmin + (xmax - xmin)*real(i - 1, dp) / (n_points - 1)
        fx(i) = func(x(i))
        dfx(i) = dfunc(x(i))
        ddfx(i) = ddfunc(x(i))
        write(1, *) x(i), fx(i), dfx(i), ddfx(i)
    end do
    write(1, *) ' '
    write(1, *) ' '

    
    ! ------------------------ Derivatives ------------------------
    call forward_derivate(n_points, x, df_for, func)
    call backward_derivate(n_points, x, df_back, func)
    call central_derivate(n_points, x, df_cen, func)
    call five_points_derivate(n_points, x, df_5p, func)

    write(1, *) '# Numeric Derivatives'
    do i = 1, n_points
        write(1, *) x(i), dfx(i), df_for(i), df_back(i), df_cen(i), df_5p(i)
    end do
    write(1, *) ' '
    write(1, *) ' '

    call second_derivate(n_points, x, ddf_num, func)
    write(1, *) '# Numeric 2nd Derivative'
    do i = 1, n_points
        write(1, *) x(i), ddfx(i), ddf_num(i)
    end do
    write(1, *) ' '
    write(1, *) ' '
    
    deallocate(x, fx, dfx, ddfx, df_for, df_back, df_cen, df_5p, ddf_num)

    call system('gnuplot -p Derivatives_fig1.gnu')
    call system('gnuplot -p Derivatives_fig2.gnu')
    call system('gnuplot -p Derivatives_fig3.gnu')
end program Derivatives

function func(x) result(fx)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: x, fx

    fx = x*exp(-x)
end function func

function dfunc(x) result(fx)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: x, fx

    fx = exp(-x) - x*exp(-x)
end function dfunc

function ddfunc(x) result(fx)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: x, fx

    fx = -2*exp(-x) + x*exp(-x)
end function ddfunc

! ------------------------ First Derivatives ------------------------
subroutine central_derivate(n_points, x, dfunc, func)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: x(n_points), dfunc(n_points), func, step
    
    step = x(2) - x(1)
    dfunc(1) = (func(x(2)) - func(x(1))) / step
    dfunc(n_points) = (func(x(n_points)) - func(x(n_points - 1))) / step
    do i = 2, n_points - 1
        dfunc(i) = (func(x(i + 1)) - func(x(i - 1))) / (2 * step)
    end do
end subroutine central_derivate

subroutine forward_derivate(n_points, x, dfunc, func)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: x(n_points), dfunc(n_points), func, step
    
    step = x(2) - x(1)
    dfunc(n_points) = (func(x(n_points)) - func(x(n_points - 1))) / step
    do i = 1, n_points - 1
        dfunc(i) = (func(x(i + 1)) - func(x(i))) / step
    end do
end subroutine forward_derivate

subroutine backward_derivate(n_points, x, dfunc, func)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: x(n_points), dfunc(n_points), func, step
    
    step = x(2) - x(1)
    dfunc(1) = (func(x(2)) - func(x(1))) / step
    do i = 2, n_points
        dfunc(i) = (func(x(i)) - func(x(i - 1))) / step
    end do
end subroutine backward_derivate

subroutine five_points_derivate(n_points, x, dfunc, func)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: x(n_points), dfunc(n_points), func, step
    
    step = x(2) - x(1)
    dfunc(1) = (func(x(2)) - func(x(1))) / step
    dfunc(2) = (func(x(3)) - func(x(1))) / (2 * step)
    dfunc(n_points) = (func(x(n_points)) - func(x(n_points - 1))) / step
    dfunc(n_points - 1) = (func(x(n_points)) - func(x(n_points - 2))) / (2 * step)
    do i = 3, n_points - 2
        dfunc(i) = (func(x(i - 2)) - func(x(i + 2)) + 8*(func(x(i + 1)) - func(x(i - 1)))) / (12 * step)
    end do
end subroutine five_points_derivate

! ------------------------ Second Derivatives ------------------------
subroutine second_derivate(n_points, x, dfunc, func)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points
    real(dp) :: x(n_points), dfunc(n_points), func, step
    
    step = x(2) - x(1)
    dfunc(1) = (func(x(3)) + func(x(1)) - 2*func(x(2))) / (step**2)
    dfunc(n_points) = (func(x(n_points)) + func(x(n_points - 2)) - 2*func(x(n_points - 1))) / (step**2)
    do i = 2, n_points - 1
        dfunc(i) = (func(x(i + 1)) + func(x(i - 1)) - 2*func(x(i))) / (step**2)
    end do
end subroutine second_derivate