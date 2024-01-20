program EDO
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: ndades, i, j, steps(4)
    real(dp) :: m, l, g, Tn, omega_n, pi, tini, tend, phi0, omega_0, step, t, &
                angular_accel, angular_accel_approx, ecin_e, ecin_e2, ecin_ab, epot_e, epot_e2, epot_ab, &
                etot_e, etot_e2, etot_ab, ecine, epoten
    real(dp), allocatable :: phi_e(:), omega_e(:), phi_e2(:), omega_e2(:), phi_ab(:), omega_ab(:), phi_a3p(:), omega_a3p(:), &
                             phi_a4p(:), omega_a4p(:), phi_am(:), omega_am(:), phi_h(:), omega_h(:)
    external :: angular_accel, angular_accel_approx

    pi = acos(-1._dp)
    m = 0.51_dp ! kg
    l = 0.45_dp ! m
    g = 3.71_dp ! m/s^2
    omega_n = sqrt(g/l)
    Tn = 2._dp * pi / omega_n
    tini = 0._dp
    tend = 6._dp * Tn

    open(1, file='EDO_res.dat')

    ! --- Small Oscilations --- 
    write(1, *) '# Simple Pendulum: Small Oscilations'
    write(1, *) '# t', 'Euler', 'Enhanced Euler', 'Adams Bashforth'

    ndades = 1300
    allocate(phi_e(ndades), omega_e(ndades), phi_e2(ndades), omega_e2(ndades), phi_ab(ndades), omega_ab(ndades))

    ! We use the degree reduced ecuation phi(t), phi'(t)=omega
    omega_0 = 0._dp
    phi0 = 0.02_dp
    step = (tend - tini) / real(ndades - 1, dp)
    call EulerMethod(ndades, step, phi0, omega_0, phi_e, omega_e, angular_accel_approx)
    call EnhancedEulerMethod(ndades, step, phi0, omega_0, phi_e2, omega_e2, angular_accel_approx)
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel_approx)
    
    do i = 1, ndades
        t = tini + (i - 1) * step
        write(1,*) t, phi_e(i), omega_e(i), phi_e2(i), omega_e2(i), phi_ab(i), omega_ab(i)
    end do
    write(1, *) ' '
    write(1, *) ' '
    deallocate(phi_e, omega_e, phi_e2, omega_e2, phi_ab, omega_ab)

    ! --- Wide Oscilations --- 
    write(1, *) '# Simple Pendulum: Wide Oscilations'
    write(1, *) '# t', 'Euler', 'Enhanced Euler', 'Adams Bashforth'

    ndades = 1800
    allocate(phi_e(ndades), omega_e(ndades), phi_e2(ndades), omega_e2(ndades), phi_ab(ndades), omega_ab(ndades), &
             phi_a3p(ndades), omega_a3p(ndades), phi_a4p(ndades), omega_a4p(ndades), phi_am(ndades), omega_am(ndades), &
             phi_h(ndades), omega_h(ndades))

    omega_0 = 0._dp
    phi0 = pi - 0.025_dp
    step = (tend - tini) / real(ndades - 1, dp)
    call EulerMethod(ndades, step, phi0, omega_0, phi_e, omega_e, angular_accel)
    call EnhancedEulerMethod(ndades, step, phi0, omega_0, phi_e2, omega_e2, angular_accel)
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel)
    call Adams_3p(ndades, step, phi0, omega_0, phi_a3p, omega_a3p, angular_accel)
    call Adams_4p(ndades, step, phi0, omega_0, phi_a4p, omega_a4p, angular_accel)
    call Adams_Moulton(ndades, step, phi0, omega_0, phi_am, omega_am, angular_accel)
    call Hamming(ndades, step, phi0, omega_0, phi_h, omega_h, angular_accel)
    
    do i = 1, ndades
        t = tini + (i - 1) * step
        write(1,*) t, phi_e(i), omega_e(i), phi_e2(i), omega_e2(i), phi_ab(i), omega_ab(i), phi_a3p(i), omega_a3p(i), &
                   phi_a4p(i), omega_a4p(i), phi_am(i), omega_am(i), phi_h(i), omega_h(i)
    end do
    write(1, *) ' '
    write(1, *) ' '
    deallocate(phi_e, omega_e, phi_e2, omega_e2, phi_ab, omega_ab, phi_a3p, omega_a3p, phi_a4p, omega_a4p, &
               phi_am, omega_am, phi_h, omega_h)

    ! --- System Energy Analysis ---
    write(1, *) '# Simple Pendulum: System Energy Analysis'
    write(1, *) '# t', 'Euler', 'Enhanced Euler', 'Adams Bashforth'
    
    ndades = 2500
    allocate(phi_e(ndades), omega_e(ndades), phi_e2(ndades), omega_e2(ndades), phi_ab(ndades), omega_ab(ndades))
    
    write(1, *) "# phi0 = pi - 0.042; phi'0 = 0"
    omega_0 = 0._dp
    phi0 = pi - 0.042_dp
    step = (tend - tini) / real(ndades - 1, dp)
    call EulerMethod(ndades, step, phi0, omega_0, phi_e, omega_e, angular_accel)
    call EnhancedEulerMethod(ndades, step, phi0, omega_0, phi_e2, omega_e2, angular_accel)
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel)
    
    do i = 1, ndades
        t = tini + (i - 1) * step

        ecin_e = Ecine(omega_e(i))
        ecin_e2 = Ecine(omega_e2(i))
        ecin_ab = Ecine(omega_ab(i))

        epot_e = Epoten(phi_e(i))
        epot_e2 = Epoten(phi_e2(i))
        epot_ab = Epoten(phi_ab(i))

        etot_e = ecin_e + epot_e
        etot_e2 = ecin_e2 + epot_e2
        etot_ab = ecin_ab + epot_ab

        write(1,*) t, epot_e, epot_e2, epot_ab, ecin_e, ecin_e2, ecin_ab, etot_e, etot_e2, etot_ab
    end do
    write(1, *) ' '
    write(1, *) ' '

    write(1, *) "# phi0 = 1; phi'0 = 0"
    phi0 = 1._dp
    call EulerMethod(ndades, step, phi0, omega_0, phi_e, omega_e, angular_accel)
    call EnhancedEulerMethod(ndades, step, phi0, omega_0, phi_e2, omega_e2, angular_accel)
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel)
    
    do i = 1, ndades
        t = tini + (i - 1) * step

        ecin_e = Ecine(omega_e(i))
        ecin_e2 = Ecine(omega_e2(i))
        ecin_ab = Ecine(omega_ab(i))

        epot_e = Epoten(phi_e(i))
        epot_e2 = Epoten(phi_e2(i))
        epot_ab = Epoten(phi_ab(i))

        etot_e = ecin_e + epot_e
        etot_e2 = ecin_e2 + epot_e2
        etot_ab = ecin_ab + epot_ab

        write(1,*) t, epot_e, epot_e2, epot_ab, ecin_e, ecin_e2, ecin_ab, etot_e, etot_e2, etot_ab
    end do
    write(1, *) ' '
    write(1, *) ' '
    deallocate(phi_e, omega_e, phi_e2, omega_e2, phi_ab, omega_ab)

    ! --- Transition ---
    write(1, *) '# Simple Pendulum: Transition'
    write(1, *) '# t', 'Euler', 'Enhanced Euler', 'Adams Bashforth'

    ndades = 2100
    allocate(phi_e(ndades), omega_e(ndades), phi_ab(ndades), omega_ab(ndades))

    phi0 = 0._dp
    tend = 7._dp * Tn
    step = (tend - tini) / real(ndades - 1, dp)

    omega_0 = 2._dp * sqrt(g/l) + 0.04_dp
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel)
    
    omega_0 = 2._dp * sqrt(g/l) - 0.04_dp
    call Adams_Bashforth(ndades, step, phi0, omega_0, phi_e, omega_e, angular_accel)
    
    do i = 1, ndades
        t = tini + (i - 1) * step
        write(1,*) t, phi_e(i), omega_e(i), phi_ab(i), omega_ab(i)
    end do
    write(1, *) ' '
    write(1, *) ' '
    deallocate(phi_e, omega_e, phi_ab, omega_ab)

    ! --- Method Convergence ---
    write(1, *) '# Simple Pendulum: Method Convergence'
    write(1, *) '# t', 'Euler', 'Enhanced Euler', 'Adams Bashforth'

    omega_0 = 0.1_dp
    phi0 = 2.1_dp
    tend = 12._dp * Tn

    steps = [300, 1000, 2200, 14500]
    do i = 1, size(steps)
        ndades = steps(i)
        step = (tend - tini) / real(ndades - 1, dp)
        write(1, *) '# ', steps(i), 'steps'
        allocate(phi_ab(ndades), omega_ab(ndades))

        call Adams_Bashforth(ndades, step, phi0, omega_0, phi_ab, omega_ab, angular_accel)
        
        do j = 1, ndades
            t = tini + (j - 1) * step
    
            ecin_ab = Ecine(omega_ab(j))
    
            epot_ab = Epoten(phi_ab(j))
    
            etot_ab = ecin_ab + epot_ab
    
            write(1,*) t, epot_ab, ecin_ab, etot_ab
        end do
        write(1, *) ' '
        write(1, *) ' '
        deallocate(phi_ab, omega_ab)
    end do

    close(1)

    call system('gnuplot -p EDO_fig1.gnu')
    call system('gnuplot -p EDO_fig2.gnu')
    call system('gnuplot -p EDO_fig3.gnu')
    call system('gnuplot -p EDO_fig4.gnu')
    call system('gnuplot -p EDO_fig5.gnu')
    call system('gnuplot -p EDO_fig6.gnu')
    call system('gnuplot -p EDO_fig7.gnu')
    call system('gnuplot -p EDO_fig8.gnu')
    call system('gnuplot -p EDO_fig8b.gnu')
end program EDO

! ---------------------- EDO Solving Methods ----------------------
subroutine EulerMethod(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step
    external :: func

    ! Euler Method: y_n = y_n-1 + h*f_n-1
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Euler method for both

    phi(1) = phi0 + step * omega_0 ! f_n = p; y_n = phi
    omega(1) = omega_0 + step * func(phi0) ! f_n = p' = -g * sin(phi) / l; y_n = p = phi'

    do i = 2, ndades
        phi(i) = phi(i - 1) + step * omega(i - 1) ! f_n = p; y_n = phi
        omega(i) = omega(i - 1) + step * func(phi(i - 1)) ! f_n = p' = -g * sin(phi) / l; y_n = p = phi'
    end do
end subroutine EulerMethod

subroutine EnhancedEulerMethod(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step
    external :: func
    
    ! Euler Method+: y_n = y_n-2 + 2*h*f_n-1
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Euler method+ for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Euler method+ with the initial state 
    phi(2) = phi0 + 2 * step * omega(1)
    omega(2) = omega_0 + 2 * step * func(phi(1))

    do i = 3, ndades
        phi(i) = phi(i - 2) + 2 * step * omega(i - 1)
        omega(i) = omega(i - 2) + 2 * step * func(phi(i - 1))
    end do
end subroutine EnhancedEulerMethod

subroutine Adams_Bashforth(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step
    external :: func

    ! Adams-Bashforth: y_n = y_n-1 + h*f_n-2 / 2 + 3*h*f_n-1 / 2
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Adams-Bashforth for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Adams-Bashforth with the initial state 
    phi(2) = phi(1) + 3 * step * omega(1) / 2._dp - step * omega_0 / 2._dp
    omega(2) = omega(1) + 3 * step * func(phi(1)) / 2._dp - step * func(phi0) / 2._dp

    do i = 3, ndades
        phi(i) = phi(i - 1) + step*(3 * omega(i - 1) - omega(i - 2)) / 2._dp
        omega(i) = omega(i - 1) + step*(3 * func(phi(i - 1)) - func(phi(i - 2))) / 2._dp
    end do
end subroutine Adams_Bashforth

subroutine Adams_3p(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step
    external :: func

    ! Adams-Bashforth: y_n = y_n-1 + h*f_n-2 / 2 + 3*h*f_n-1 / 2
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Adams-Bashforth for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Euler method+ with the initial state 
    phi(2) = phi0 + 2 * step * omega(1)
    omega(2) = omega_0 + 2 * step * func(phi(1))

    ! Lastly we apply the Adams_3p method with the initial state 
    phi(3) = phi(2) + step * (23 * omega(2) - 16 * omega(1) + 5 * omega_0) / 12._dp
    omega(3) = omega(2) + step * (23 * func(phi(2)) - 16 * func(phi(1)) + 5 * func(phi0)) / 12._dp

    do i = 4, ndades
        phi(i) = phi(i - 1) + step * (23 * omega(i - 1) - 16 * omega(i - 2) + 5 * omega(i - 3)) / 12._dp
        omega(i) = omega(i - 1) + step * (23 * func(phi(i - 1)) - 16 * func(phi(i - 2)) + 5 * func(phi(i - 3))) / 12._dp
    end do
end subroutine Adams_3p

subroutine Adams_4p(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step
    external :: func

    ! Adams-Bashforth: y_n = y_n-1 + h*f_n-2 / 2 + 3*h*f_n-1 / 2
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Adams-Bashforth for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Euler method+ with the initial state 
    phi(2) = phi0 + 2 * step * omega(1)
    omega(2) = omega_0 + 2 * step * func(phi(1))

    ! Then we apply the Adams_3p method with the initial state 
    phi(3) = phi(2) + step * (23 * omega(2) - 16 * omega(1) + 5 * omega_0) / 12._dp
    omega(3) = omega(2) + step * (23 * func(phi(2)) - 16 * func(phi(1)) + 5 * func(phi0)) / 12._dp

    ! Lastly we apply the Adams_4p method with the initial state 
    phi(4) = phi(3) + step * (55 * omega(3) - 59 * omega(2) + 37 * omega(1) - 9 * omega_0) / 24._dp
    omega(4) = omega(3) + step * (55 * func(phi(3)) - 59 * func(phi(2)) + 37 * func(phi(1)) - 9 * func(phi0)) / 24._dp

    do i = 5, ndades
        phi(i) = phi(i - 1) + step * (55 * omega(i - 1) - 59 * omega(i - 2) + 37 * omega(i - 3) - 9 * omega(i - 4)) / 24._dp
        omega(i) = omega(i - 1) + step * (55 * func(phi(i - 1)) - 59 * func(phi(i - 2)) + 37 * func(phi(i - 3)) &
                   - 9 * func(phi(i - 4))) / 24._dp
    end do
end subroutine Adams_4p

subroutine Adams_Moulton(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step, phi_p, omega_p
    external :: func

    ! Adams-Bashforth: y_n = y_n-1 + h*f_n-2 / 2 + 3*h*f_n-1 / 2
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Adams-Bashforth for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Euler method+ with the initial state 
    phi(2) = phi0 + 2 * step * omega(1)
    omega(2) = omega_0 + 2 * step * func(phi(1))

    ! Then we apply the Adams_3p method with the initial state 
    phi(3) = phi(2) + step * (23 * omega(2) - 16 * omega(1) + 5 * omega_0) / 12._dp
    omega(3) = omega(2) + step * (23 * func(phi(2)) - 16 * func(phi(1)) + 5 * func(phi0)) / 12._dp

    ! Lastly we apply the Adams_Moulton method with the initial state 
    phi_p = phi(3) + step * (55 * omega(3) - 59 * omega(2) + 37 * omega(1) - 9 * omega_0) / 24._dp
    omega_p = omega(3) + step * (55 * func(phi(3)) - 59 * func(phi(2)) + 37 * func(phi(1)) - 9 * func(phi0)) / 24._dp

    phi(4) = phi(3) + step * (9 * omega_p + 19 * omega(3) - 5 * omega(2) + omega(1)) / 24._dp
    omega(4) = omega(3) + step * (9 * func(phi_p) + 19 * func(phi(3)) - 5 * func(phi(2)) + func(phi(1))) / 24._dp

    do i = 5, ndades
        phi_p = phi(i - 1) + step * (55 * omega(i - 1) - 59 * omega(i - 2) + 37 * omega(i - 3) - 9 * omega(i - 4)) / 24._dp
        omega_p = omega(i - 1) + step * (55 * func(phi(i - 1)) - 59 * func(phi(i - 2)) + 37 * func(phi(i - 3)) &
                  - 9 * func(phi(i - 4))) / 24._dp

        phi(i) = phi(i - 1) + step * (9 * omega_p + 19 * omega(i - 1) - 5 * omega(i - 2) + omega(1)) / 24._dp
        omega(i) = omega(i - 1) + step * (9 * func(phi_p) + 19 * func(phi(i - 1)) - 5 * func(phi(i - 2)) &
                   + func(phi(i - 3))) / 24._dp
    end do
end subroutine Adams_Moulton

subroutine Hamming(ndades, step, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, ndades
    real(dp):: phi0, omega_0, phi(ndades), omega(ndades), func, step, phi_p, omega_p, phi_m, omega_m, phi_c, omega_c, &
               phi_p_prev, omega_p_prev, phi_c_prev, omega_c_prev
    external :: func

    ! Adams-Bashforth: y_n = y_n-1 + h*f_n-2 / 2 + 3*h*f_n-1 / 2
    ! We got two eqs -> phi = omega; alpha = -g·sin(phi)/l (see angular_accel for deduction)
    ! Apply Adams-Bashforth for both

    ! For the first point we apply simple Euler method
    phi(1) = phi0 + step * omega_0
    omega(1) = omega_0 + step * func(phi0)

    ! Then we apply the Euler method+ with the initial state 
    phi(2) = phi0 + 2 * step * omega(1)
    omega(2) = omega_0 + 2 * step * func(phi(1))

    ! Then we apply the Adams_3p method with the initial state 
    phi(3) = phi(2) + step * (23 * omega(2) - 16 * omega(1) + 5 * omega_0) / 12._dp
    omega(3) = omega(2) + step * (23 * func(phi(2)) - 16 * func(phi(1)) + 5 * func(phi0)) / 12._dp

    ! Lastly we apply the Adams_Moulton method with the initial state and y_p_prev = y_c_prev
    phi_p = phi0 + 4 * step * (2 * omega(3) - omega(2) + 2 * omega(1)) / 3._dp
    omega_p = omega_0 + 4 * step * (2 * func(phi(3)) - func(phi(2)) + 2 * func(phi(1))) / 3._dp

    phi_m = phi_p
    omega_m = omega_p

    phi_c = (9 * phi(3) - phi(1) + 3 * step * (omega_m + 2 * omega(3) - omega(2))) / 8._dp
    omega_c = (9 * omega(3) - omega(1) + 3 * step * (func(phi_m) + 2 * func(phi(3)) - func(phi(2)))) / 8._dp

    phi(4) = phi_c
    omega(4) = omega_c

    do i = 5, ndades
        phi_p = phi(i - 4) + 4 * step * (2 * omega(i - 1) - omega(i - 2) + 2 * omega(i - 3)) / 3._dp
        omega_p = omega(i - 4) + 4 * step * (2 * func(phi(i - 1)) - func(phi(i - 2)) + 2 * func(phi(i - 3))) / 3._dp

        phi_m = phi_p - 112 * step * (phi_p_prev - phi_c_prev) / 121._dp
        omega_m = omega_p - 112 * step * (omega_p_prev - omega_c_prev) / 121._dp

        phi_c = (9 * phi(i - 1) - phi(i - 3) + 3 * step * (omega_m + 2 * omega(i - 1) - omega(i - 2))) / 8._dp
        omega_c = (9 * omega(i - 1) - omega(i - 3) + 3 * step * (func(phi_m) + 2 * func(phi(i - 1)) - func(phi(i - 2)))) / 8._dp

        phi(i) = phi_c + 9 * (phi_p - phi_c) / 121._dp
        omega(i) = omega_c + 9 * (omega_p - omega_c) / 121._dp
        
        phi_p_prev = phi_p
        omega_p_prev = omega_p
        phi_c_prev = phi_c
        omega_c_prev = omega_c
    end do
end subroutine Hamming

! ---------------------- EDO Presentation ----------------------
function angular_accel(phi) result(alpha)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp), intent(in) :: phi
    real(dp) :: alpha, g, l

    l = 0.45_dp ! cm
    g = 3.71_dp ! m/s^2

    ! Original: l·phi''(t) = -g·sin(phi(t))
    ! Degree reduced: omega = phi'; omega' = phi'' -> l·omega'(t) = -g·sin(phi(t)) 
    ! omega'(t) = alpha = -g·sin(phi(t)) / l

    alpha = -g*sin(phi) / l
end function angular_accel

function angular_accel_approx(phi) result(alpha)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp), intent(in) :: phi
    real(dp) :: alpha, g, l

    l = 0.45_dp ! cm
    g = 3.71_dp ! m/s^2

    alpha = -g*phi / l
end function angular_accel_approx

function Ecine(omega) result(ecin)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp), intent(in) :: omega
    real(dp) :: ecin, g, l, pi, m

    pi = acos(-1._dp)
    l = 0.45_dp ! cm
    g = 3.71_dp ! m/s^2
    m = 0.51_dp ! kg

    ecin = m*((omega*l)**2) / 2._dp
end function Ecine

function Epoten(phi) result(epot)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp), intent(in) :: phi
    real(dp) :: epot, g, l, pi, m

    pi = acos(-1._dp)
    l = 0.45_dp ! cm
    g = 3.71_dp ! m/s^2
    m = 0.51_dp ! kg

    epot = -m*g*l*cos(phi)
end function Epoten