program gener2024
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/gamma, x0, pi, lambda, omega

    integer :: i, j, N, n_points
    real(dp) :: gamma, x0, I1, pi, xini, xend, integral, func, std, step, lambda, omega, &
                phi0, zeta_n(2), zeta0, dzeta, dphi, energy
    real(dp), allocatable :: x_out(:), t(:), zeta(:), phi(:)
    external :: I1, func, dzeta, dphi, energy

    ! ---------------- Problem 1 ----------------
    ! P1 a)
    open(1, file='Exa-jan-24-res.dat')
    n_points = 1000
    step = 10._dp / real(n_points - 1, dp)
    allocate(t(n_points), zeta(n_points), phi(n_points))

    ! --- Generate t ---
    do i = 1, n_points
        t(i) = step*(i-1)
    end do

    lambda = 2.5_dp
    omega = 1._dp
    phi0 = 0._dp
    zeta_n = [0.3_dp, 0.9_dp]

    do i = 1, size(zeta_n)
        ! --- Differential equation solving via RK4 ---
        zeta0 = zeta_n(i)
        write(1, *) '# Initial Conditions zeta0=', zeta0, '; phi0=', phi0
        call RungeKutta4order(n_points, step, t, zeta0, phi0, zeta, phi, dzeta, dphi)
        do j = 1, n_points
            write(1, '(3f20.14)') t(j), zeta(j), phi(j)
        end do
        write(1, *) ' '
        write(1, *) ' '

        ! --- Energy comparison between t0 and tEnd ---
        write(1, *) '# Energy at the t0 vs tEnd'
        write(1, '(2f20.14)') t(1), energy(zeta(1), phi(1))
        write(1, '(2f20.14)') t(n_points), energy(zeta(n_points), phi(n_points))
        print *, 'Energy at t =', t(1), '; E =', energy(zeta(1), phi(1))
        print *, 'Energy at t =', t(n_points), '; E =', energy(zeta(n_points), phi(n_points))
        write(1, *) ' '
        write(1, *) ' '
    end do
    
    deallocate(t, zeta, phi)
    
    ! ---------------- Problem 2 ----------------
    pi = acos(-1._dp)
    x0 = 1._dp
    gamma = 0.5_dp

    ! P2 a)
    xini = 0._dp
    xend = pi

    ! --- Use of Simpsons rule to integrate I1 ---
    call simpson(xini, xend, 15, integral, I1)
    print *, 'Integral I1 = ', integral
    
    ! P2 b)
    ! --- Generation of random samples --- 
    allocate(x_out(100000))
    call box_muller(100000, x_out, 0._dp, 1._dp) 
    ! Box-Muller is not the correct way. Probably a hybrid version between AR and inversion would be best.
    ! For example generating sin(x)**2 via inversion and using it to envelope I1 to apply AR.
    
    write(1, *) '# Montecarlo integral'
    do i = 1, 10
        N = i*10000
        
        call weightedmontecarlo1D(N, x_out, integral, std, I1)
        write(1, *) real(N*10000, dp), integral, std
        print *, 'Integral MC I1 = ', integral
    end do
    write(1, *) ' '
    write(1, *) ' '

    deallocate(x_out)

    close(1)

    call system('gnuplot -p Exa-jan-24-fig1.gnu')
    call system('gnuplot -p Exa-jan-24-fig2.gnu')
end program gener2024

! ---------------- Problem 1 ----------------
function dzeta(t, zeta, phi) result(fzeta)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/gamma, x0, pi, lambda, omega

    real(dp) :: x, fx, gamma, x0, pi, lambda, omega, t, zeta, phi, fzeta

    fzeta = -omega*sqrt(1._dp-(zeta**2))*sin(phi)
end function dzeta

function dphi(t, zeta, phi) result(fphi)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/gamma, x0, pi, lambda, omega

    real(dp) :: gamma, x0, pi, lambda, omega, t, zeta, phi, fphi

    fphi = lambda*zeta + omega*zeta*cos(phi)/sqrt(1._dp-(zeta**2))
end function dphi

function energy(zeta, phi) result(E)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/gamma, x0, pi, lambda, omega

    real(dp) :: gamma, x0, pi, lambda, omega, zeta, phi, E

    E = lambda * (zeta**2)/2._dp - sqrt(1._dp - zeta**2)*cos(phi)*omega
end function energy

subroutine RungeKutta4order(n_points, step, t, zeta0, phi0, zeta, phi, fzeta, fphi)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i, j
    real(dp):: step, k1(2), k2(2), k3(2), k4(2), phi0, zeta0, &
               phi(n_points), zeta(n_points), fzeta, fphi, t(n_points) ! x for V = V(x)
    external :: fzeta, fphi

    ! We got two eqs -> zeta' = fzeta; phi' = fphi
    
    ! --- RK4 ---
    ! k1 = f(x0, y0)
    ! k2 = f(x0 + h/2, y0 + h·k1/2)
    ! k3 = f(x0 + h/2, y0 + h·k2/2)
    ! k4 = f(x0 + h, y0 + h·k3)
    ! y1 = y0 + h(k1 + 2·k2 + 2·k3 + k4)/6
    ! ------------

    k1(1) = fzeta(t(1), zeta0, phi0)
    k1(2) = fphi(t(1), zeta0, phi0)

    k2(1) = fzeta(t(1) + step/2._dp, zeta0 + step*k1(1)/2._dp, phi0 + step*k1(2)/2._dp)
    k2(2) = fphi(t(1) + step/2._dp, zeta0 + step*k1(1)/2._dp, phi0 + step*k1(2)/2._dp)

    k3(1) = fzeta(t(1) + step/2._dp, zeta0 + step*k2(1)/2._dp, phi0 + step*k2(2)/2._dp)
    k3(2) = fphi(t(1) + step/2._dp, zeta0 + step*k2(1)/2._dp, phi0 + step*k2(2)/2._dp)

    k4(1) = fzeta(t(1) + step, zeta0 + step*k3(1), phi0 + step*k3(2))
    k4(2) = fphi(t(1) + step, zeta0 + step*k3(1), phi0 + step*k3(2))

    zeta(1) = zeta0 + step*(k1(1) + 2._dp*k2(1) + 2._dp*k3(1) + k4(1))/6._dp
    phi(1) = phi0 + step*(k1(2) + 2._dp*k2(2) + 2._dp*k3(2) + k4(2))/6._dp

    do i = 2, n_points
        k1(1) = fzeta(t(i), zeta0, phi0)
        k1(2) = fphi(t(i), zeta0, phi0)

        k2(1) = fzeta(t(i) + step/2._dp, zeta(i-1) + step*k1(1)/2._dp, phi(i-1) + step*k1(2)/2._dp)
        k2(2) = fphi(t(i) + step/2._dp, zeta(i-1) + step*k1(1)/2._dp, phi(i-1) + step*k1(2)/2._dp)

        k3(1) = fzeta(t(i) + step/2._dp, zeta(i-1) + step*k2(1)/2._dp, phi(i-1) + step*k2(2)/2._dp)
        k3(2) = fphi(t(i) + step/2._dp, zeta(i-1) + step*k2(1)/2._dp, phi(i-1) + step*k2(2)/2._dp)

        k4(1) = fzeta(t(i) + step, zeta(i-1) + step*k3(1), phi(i-1) + step*k3(2))
        k4(2) = fphi(t(i) + step, zeta(i-1) + step*k3(1), phi(i-1) + step*k3(2))

        zeta(i) = zeta(i-1) + step*(k1(1) + 2._dp*k2(1) + 2._dp*k3(1) + k4(1))/6._dp
        phi(i) = phi(i-1) + step*(k1(2) + 2._dp*k2(2) + 2._dp*k3(2) + k4(2))/6._dp
    end do
end subroutine RungeKutta4order

! ---------------- Problem 2 ----------------
function I1(x) result(fx)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/gamma, x0, pi, lambda, omega

    real(dp) :: x, fx, gamma, x0, pi, lambda, omega

    fx = (sin(x)**2)/((x-x0)**2 + (gamma)**2)
end function I1

subroutine simpson(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, a, fa, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = funcio(xend)

    do i = 1, n_iter - 1
        a = real(i, dp) * step + xini
        if (mod(i, 2) .eq. 0) then
            fa = funcio(a)
            AS = AS + 2 * fa
        else
            fa = funcio(a)
            AS = AS + 4 * fa
        end if
    end do
    AS = step * AS / 3._dp
end subroutine simpson

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
    std = sqrt(std / real(N, dp) - integral**2) / sqrt(real(N, dp))
end subroutine weightedmontecarlo1D

subroutine box_muller(n_points, x_out, mu, sigma)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i
    real(dp) :: x_out(n_points), x1, x2, r, phi, pi, sigma, mu

    pi = acos(-1._dp)

    do i = 1, n_points, 2
        call random_number(x1)
        call random_number(x2)
        x1 = mu + x1*sigma
        x2 = mu + x2*sigma
        r = sqrt(-2._dp*log(x1))
        phi = 2._dp * pi * x2
        x_out(i) = r * cos(phi)
        x_out(i+1) = r * sin(phi)
    end do
end subroutine box_muller