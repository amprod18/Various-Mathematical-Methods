program Montecarlo_Integration
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: N, i, ndades, iseed
    real(dp) :: x_min, x_max, std1, std2, std3_ar, std3_m, integral1, integral2, integral3_ar, &
                integral3_m, cotasup, distr1, distr2, I1, I2, L, I1_real, x_distr(1000000), pi
    external :: I1, I2, distr1, distr2

    pi = acos(-1._dp)
    iseed = 20347972
    call srand(iseed)

    open(1, file='Montecarlo_Integration_res.dat')
    write(1, *) "#Ex1"
    write(1, *) "#N",  "#I1", "#Std1", "#Error1"

    ! Integral de Montecarlo 1D
    ! a) Integral 1

    x_min = 0._dp
    x_max = 2._dp * pi
    I1_real = (pi**2) * (2._dp * pi**2 + 3._dp/2._dp)
    do i = 1, 100
        N = i*10000
        call montecarlo1D(N, x_min, x_max, integral1, std1, I1)
        write(1, *) N, integral1, std1, abs(integral1 - I1_real)
    end do
    write(1, *) ' '
    write(1, *) ' '

    ! b) 
    write(1, *) "#N",  "#I2", "#Std2"
    ndades = 1000000
    L = 50._dp ! e-6
    x_min = 0._dp
    x_max = 2._dp*L
    cotasup = 0.2_dp

    call accept_reject(ndades, x_distr, x_min, x_max, cotasup, distr1)

    ! c) Integral 2 
    do i = 1, 100
        N = i*10000
        call weightedmontecarlo1D(N, x_distr, integral2, std2, I2)
        write(1, *) N, integral2, std2
    end do
    write(1, *) ' '
    write(1, *) ' '

    ! Fermions
    ! Integral 3
    write(1, *) "#Ex2"
    write(1, *) "#N",  "#I2", "#Std3"
    
    do i = 1, 30
        N = i*10000
        call ar_integration(N, integral3_ar, x_min, x_max, cotasup, std3_ar, I2)
        call weightedmontecarlo1D(N, x_distr, integral3_m, std3_m, I2)
        write(1, *) N, integral3_ar, std3_ar, integral3_m, std3_m
    end do
    write(1, *) ' '
    write(1, *) ' '

    close(1)

    call system('gnuplot -p Montecarlo_Integration_fig1.gnu')
    call system('gnuplot -p Montecarlo_Integration_fig2.gnu')
    call system('gnuplot -p Montecarlo_Integration_fig3.gnu')
end program Montecarlo_Integration

! ---------------------- Montecarlo Integration ----------------------
subroutine montecarlo1D(N, x_min, x_max, integral, std, func)
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
end subroutine montecarlo1D

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

subroutine ar_integration(N, area, xlow, xhigh, cotasup, std, funcio)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: N, i, count
    real(dp) :: xlow, xhigh, cotasup, funcio, px, area, x1, x2, p, x, std
    external :: funcio

    ! Area calculation via accept-reject
    do i = 1, N
        call random_number(x1)
        call random_number(x2)
        p = cotasup*x2
        x = (xhigh - xlow)*x1 + xlow
        px = funcio(x)
        if (px >= p) then
            count = count + 1
        end if
    end do
    area = cotasup * (xhigh - xlow) * real(count, dp)/real(N, dp)
    std = sqrt((area ** 2 - area / (cotasup * (xhigh - xlow))) / real(N, dp))
end subroutine ar_integration

function I1(x) result(y)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: pi, x, y

    pi = acos(-1._dp)

    y = (x**3) * cos(x)**2
end function I1

function I2(x) result(y)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: pi, x, y, distr1, distr2

    pi = acos(-1._dp)

    y = distr1(x) * distr2(x)
end function I2

function distr1(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: pi, x, px, L
    pi = acos(-1._dp)
    L = 50._dp ! e-6

    px = (sin(pi*(x - 2._dp * L) / (2._dp * L)))**2 / L
end function distr1

function distr2(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: pi, x, px, L
    pi = acos(-1._dp)
    L = 50._dp ! e-6

    px = (sin(pi*(x - 2._dp * L) / L))**2 / L
end function distr2