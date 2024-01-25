program gener2017
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, j, n_points, n_iter, N, n_caixes
    real(dp) :: step, theta, ftheta, f, xini, xend, preci, root_n(6), integral(10), step_n(10), mc_integral, std, xlow, xhigh, &
                boxsize, x, px, dist, I2, integral2, int_real, cotasup, dist2
    real(dp), allocatable :: f_int(:), x_values(:), xhis(:), vhis(:), errhis(:)
    external :: I2, dist2

    open(1, file='res1a.dat')

    n_points = 1000
    xini = -5._dp
    xend = 9._dp
    step = (xend - xini)/real(n_points - 1, dp)
    
    write(1, *) '# theta, f(theta)'
    do i = 1, n_points
        theta = xini + step*(i - 1)
        ftheta = f(theta)
        write(1, '(2f20.12)') theta, ftheta
    end do
    write(1, *) ' '
    write(1, *) ' '
    
    preci = 0.00000001_dp
    call subE1(f, preci, n_iter, root_n)
    
    write(1, *) '# Roots of f(theta)'
    do i = 1, 6
        write(1, '(f20.10)') root_n(i)
    end do
    write(1, *) ' '
    write(1, *) ' '

    write(1, *) '# Step, Int/Int_N'
    do i = 2, 20, 2
        step_n(i/2) = (xend - xini)/(2**i - 1)
        allocate(f_int(2**i))
        do j = 1, 2**i
            f_int(j) = f(5._dp + step_n(i/2)*(j - 1))
        end do
        call subE2(2**i, step_n(i/2), f_int, integral(i/2))
        deallocate(f_int)
    end do
    do i = 1, 10
        write(1, '(2f20.14)') step_n(i), abs(integral(i))
    end do
    write(1, *) ' '
    write(1, *) ' '

    open(2, file='mcE3.dat')
    
    do i = 1, 100
        N = 100*i
        call subE3(N, xini, xend, mc_integral, std, f)
        write(2, '(3f20.14)') real(N, dp), mc_integral, std
    end do
    
    N = 100000
    n_caixes = 20
    xlow = 0._dp
    xhigh = 3.5_dp
    allocate(x_values(N))
    call subE4(N, x_values)

    write(1, *) '# Histogram'
    allocate(xhis(n_caixes), vhis(n_caixes), errhis(n_caixes))
    call histograma(N, x_values, xlow, xhigh, n_caixes, xhis, vhis, errhis, boxsize)
    do i = 1, n_caixes
        x = xlow + (xhigh - xlow)*(i - 1)/real(n_caixes - 1, dp)
        px = dist(x)
        write (1, '(6f20.14)') xhis(i), vhis(i), boxsize, errhis(i), x, px
    end do
    write (1, *) ' '
    write (1, *) ' '
    deallocate(xhis, vhis, errhis, x_values)

    write(1, *) '# Integral eveloped'

    allocate(x_values(100000))
    call accept_reject(N, x_values, xlow, xhigh, cotasup, dist2)

    int_real = 27._dp/145._dp
    cotasup = 1._dp
    xhigh = 1000000000._dp
    do i = 1, 100
        N = i*1000
        call weightedmontecarlo1D(N, x_values, integral2, std, I2)
        write(1, '(4f20.14)') real(N, dp), integral2, std, abs(integral2 - int_real)
    end do
    write(1, *) ' '
    write(1, *) ' '
    deallocate(x_values)

    close(1)
    close(2)

    call system('gnuplot -p figE1.gnu')
    call system('gnuplot -p figE2.gnu')
    call system('gnuplot -p figE3.gnu')
    call system('gnuplot -p figE4.gnu')
end program gener2017

function f(x) result(fx)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: x, fx

    fx = 3._dp - abs(x - 4*cos(x))
end function f

function dist(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: x, px

    px = 3*exp(-3*x)
end function dist

function dist2(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: x, px

    px = (cos(x))**2
end function dist2

function I2(x) result(fx)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: x, fx, dist2

    fx = exp(-5*x) * dist2
end function I2

subroutine subE1(fun, preci, n_iter, root_n)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_iter, max_iter, i, j
    real(dp) :: preci, A, B, C, fa, fb, fc, fun, root_n(6), values(12)
    max_iter = 10000

    values = [-5._dp, -4._dp, -2._dp, 0._dp, 0._dp, 1._dp, 1._dp, 2._dp, 5._dp, 6._dp, 6._dp, 7._dp]

    do i = 1, 6
        A = values(2*i - 1)
        B = values(2*i)
        fa = fun(A)
        fb = fun(B)
        n_iter = ceiling(log((B - A) / preci) / log(2._dp))

        do j = 1, n_iter
            C = (A + B) / 2._dp
            fc = fun(C)
            if (fa * fb > 0) then
                print *, "Bisection method fails in the interval [", A, ", ", B, "]"
                exit
            end if
            if (fc .EQ. 0.) then
                root_n(i) = C
                exit
            end if
            if (fa * fc < 0) then
                B = C
                fb = fc
            else
                A = C
                fa = fc
            end if
            if (abs(B - A) < preci) then
                n_iter = i
                root_n(i) = C
                exit
            end if
        end do
        root_n(i) = C
    end do
end subroutine subE1

subroutine subE2(n_points, step, fx, integral)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_points, i
    real(dp) :: integral, fx(n_points), step

    integral = fx(1) + fx(n_points)

    do i = 1, n_points - 1
        if (mod(i, 2) .eq. 0) then
            integral = integral + 2 * fx(i)
        else
            integral = integral + 4 * fx(i)
        end if
    end do
    integral = step * integral / 3._dp
end subroutine subE2

subroutine subE3(N, x_min, x_max, integral, std, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: N, i
    real(dp):: x_min, x_max, integral, std, x, y, func

    integral = 0.0_dp
    std = 0.0_dp

    do i = 1, N
        call random_number(x)
        x = x_min + x * (x_max - x_min)
        y = (x_max - x_min) * func(x)
        integral = integral + y
        std = std + y * y
    end do

    integral = integral / real(N, dp)
    std = sqrt(std / real(N, dp) - integral ** 2) / sqrt(real(N, dp))
end subroutine subE3

subroutine subE4(N, x_values)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: N, i
    real(dp):: x, x_values(N)

    do i = 1, N
        call random_number(x)
        x_values(i) = -log(1 - x) / 3._dp
    end do
end subroutine subE4

subroutine histograma(ndat, xdat, xa, xb, nbox, xhis, vhis, errhis, boxsize)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: ndat, nbox, ierr, i, ibox, icount
    real(dp) :: xa, xb, boxsize
    real(dp) :: xdat(ndat), xhis(nbox), vhis(nbox), errhis(nbox)

    if (xa .ge. xb) then
        ierr = 1
        return
    end if

    boxsize = (xb - xa) / nbox

    icount = 0

    do i = 1, nbox
        vhis(i) = 0
        errhis(i) = 0
    end do

    do i = 1, ndat
        ! data between xa and xb
        if (xdat(i) .ge. xa .and. xdat(i).le.xb) then
            ibox = int((xdat(i) - xa) / boxsize) + 1 
            if (ibox .eq. nbox + 1) ibox = nbox ! puts xb into the last box
            vhis(ibox) = vhis(ibox) + 1
            icount = icount + 1
        end if
    end do
    
    if ( icount .eq. 0 ) then
        ierr = 2
        return
    end if

    ierr = 0
    print *, "Accepted: ", icount, "out of:", ndat

    do i = 1, nbox
        xhis(i) = xa + boxsize/2._dp + (i - 1)*boxsize ! central value
        errhis(i) = sqrt(vhis(i)/icount * (1._dp - vhis(i)/icount))/boxsize / sqrt(icount * 1._dp) ! errorbar
        vhis(i) = vhis(i) / icount / boxsize ! normalized value
    end do
end subroutine histograma

subroutine accept_reject(n_points, x_out, xlow, xhigh, cotasup, funcio)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i, j
    real(dp) :: xlow, xhigh, cotasup, funcio, px, x_out(n_points), x1, x2, p, x
    external :: funcio

    i = 0
    j = 1
    ! Acceptació-Rebuig
    do while(i < n_points)
        j = j + 1
        call random_number(x1)
        call random_number(x2)
        p = cotasup*x2
        x = (xhigh - xlow)*x1 + xlow
        px = funcio(x)
        if (px >= p) then
            x_out(i) = x
            i = i + 1
        end if
    end do
end subroutine accept_reject

subroutine weightedmontecarlo1D(N, x, integral, std, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: N, i
    real(dp):: integral, std, y, func, x(N)

    integral = 0.0_dp
    std = 0.0_dp

    do i = 1, N
        y = func(x(i))
        integral = integral + y
        std = std + y * y
    end do

    integral = integral / real(N, dp)
    std = sqrt(std / real(N, dp) - integral ** 2) / sqrt(real(N, dp))
end subroutine weightedmontecarlo1D