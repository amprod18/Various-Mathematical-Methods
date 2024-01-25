program EDO2
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/ E, fact, beta

    integer :: i, j, k, nequs, n_points, upper_limit, lower_limit
    real(dp) :: E, fact, xini, xend, step, omega_0, E1, E2, pi, L, delta, beta, betan(3), schrodinger_eq, schrodinger_eq_modified, &
                lowest_E
    real(dp), allocatable :: x(:), xp(:), phi0(:), phi(:, :), omega(:), phi2(:), En(:), norms(:)
    external :: RungeKutta2order, RungeKutta4order, Cash_Karp, schrodinger_eq, schrodinger_eq_modified

    pi = acos(-1._dp)
    nequs = 2
    L = 14._dp
    delta = 0.4_dp
    xini = -L/2._dp
    xend = L/2._dp
    fact = 3.80995_dp

    open(1, file='EDO2_res.dat')
    ! There's an error; we have n_points for phi (i.e.) + phi0; should be n_points - 1 for phi + phi0; same for omega
    ! ------------------------- Solving Schrodinger Eq -------------------------
    write(1, *) '# RK4 1D, Monoparticle, time independent Schrödinger Eq'
    write(1, *) '# Energy Eigenvalues and Eigenvectors'

    ! --- Testing some 'E' values ---
    n_points = 400
    step = L / real(n_points - 1, dp)
    allocate(phi0(nequs - 1), norms(nequs - 1))
    do i = 1, nequs - 1
        phi0(i) = 0._dp
    end do
    omega_0 = 0.000002_dp
    allocate(x(n_points), phi(n_points, nequs - 1), omega(n_points))

    do i = 1, n_points
        x(i) = step * real(i - 1) + xini
    end do
    
    print *, 'Final x Point:', x(n_points)/L

    allocate(En(4))
    En = [-31._dp, -30._dp, -14._dp, -13._dp]
    do i = 1, size(En)
        write(1, *)'# Not Eigenvalues'
        E = En(i)
        
        call RungeKutta4order(n_points, nequs, step, x, phi0, omega_0, phi, omega, schrodinger_eq)
        
        do j = 1, n_points
            if (x(j) .gt. L/14._dp) then
                exit
            else
                write(1,*) x(j), phi(j, :), omega(j) 
            end if
        end do
        write(1,*)' '
        write(1,*)' '
    end do
    deallocate(En)

    ! --- Finding 'E' eigenvalues ---
    allocate(En(6))
    En = [-31._dp, -30._dp, -14._dp, -13._dp, -4._dp, -3.5_dp]

    allocate(phi2(n_points))

    do i = 1, size(En), 2
        E1 = En(i)
        E2 = En(i + 1)
        call tir(n_points, nequs, x, E1, E2, phi0, omega_0, phi, omega, schrodinger_eq, RungeKutta4order) 
        print *, 'Eigenvalue Found:', E2
        E = E2 ! E2 is the resulting eigenvalue
        if (i .eq. 1) then
            lowest_E = E2
        end if

        call RungeKutta4order(n_points, nequs, step, x, phi0, omega_0, phi, omega, schrodinger_eq)
        
        do j = 1, nequs - 1
            do k = 1, n_points
                phi2(k) = phi(k, j)**2
            end do
            call simpson(n_points, step, phi2, norms(j))
            print *, norms(j)
            phi(:, j) = phi(:, j)/sqrt(norms(j))
        end do

        write(1, *) '# Eigenvectors'
        do j = 1, n_points
            write(1,*) x(j), phi(j, :), omega(j) 
        end do
        write(1,*)' '
        write(1,*)' '
    end do

    ! ------------------------- Solving Modified Schrodinger Eq -------------------------
    allocate(xp(n_points))
    betan = [0._dp, 5._dp, 15._dp]
    E = lowest_E
    
    lower_limit = ceiling((n_points - 1)*(-delta - xini)/(xend - xini) + 1)
    upper_limit = floor((n_points - 1)*(delta - xini)/(xend - xini) + 1)
    print *, lower_limit, upper_limit
    do i = 1, size(betan)
        beta = betan(i)
        phi = 0

        call RungeKutta4order(n_points, nequs, step, x, phi0, omega_0, phi, omega, schrodinger_eq_modified)
        
        do j = 1, nequs - 1
            do k = 1, n_points
                phi2(k) = phi(k, j)**2
            end do
            call simpson(n_points, step, phi2, norms(j))
            phi(:, j) = phi(:, j)/sqrt(norms(j))
        end do

        write(1, *) '# Eigenvectors Beta'
        do j = 1, n_points
            write(1,*) x(j), phi(j, :), omega(j) 
        end do
        write(1,*) ' '
        write(1,*) ' '

        do j = 1, nequs - 1
            do k = 1, n_points
                phi2(k) = phi(k, j)**2
            end do
            call simpson(upper_limit - lower_limit, step, phi2(lower_limit:upper_limit), norms(j))
        end do
        write(1, *) '# Probability of -delta <= x <= delta'
        write(1, *) beta, norms(1)
        write(1,*) ' '
        write(1,*) ' '
        print *, 'Beta =', beta, '; P(-delta <= x <= delta):', norms(1)
    end do

    deallocate(phi0, phi, phi2, omega, norms, x, xp, En)
    
    close(1)

    call system('gnuplot -p EDO2_fig1.gnu')
    call system('gnuplot -p EDO2_fig2.gnu')
    call system('gnuplot -p EDO2_fig3.gnu')
    call system('gnuplot -p EDO2_fig4.gnu')
    
end program EDO2

! ---------------------- EDO Solving Methods ----------------------
subroutine update_ki(nequs, x, y, phi0, omega_0, ki, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, nequs
    real(dp) :: x, y, func, ki(nequs), phi0(nequs), omega_0

    do i = 1, nequs - 2
        ki(i) = phi0(i + 1)
    end do
    ki(nequs - 1) = omega_0
    ki(nequs) = func(x, y)
end subroutine update_ki

subroutine RungeKutta2order(n_points, nequs, step, x, phi0, omega_0, phi, omega, func) ! Works badly for this case
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, nequs, i, j
    real(dp):: step, k1(nequs), k2(nequs), phi0(nequs - 1), omega_0, q, a1, a2, &
               phi(n_points, nequs - 1), omega(n_points), func, x(n_points) ! x for V = V(x)
    external :: func, update_ki

    ! We got two eqs -> phi' = omega; alpha = 2·phi(x)·[V(x) - E] (see schrodinger_eq for deduction)
    
    ! --- RK2 ---
    ! k1 = f(x0, y0)
    ! k2 = f(x0 + q·h, y0 + q·h·k1)
    ! y1 = y0 + h(a1·k1 + a2·k2)

    ! a1 = 1 - a2
    ! q = 1/2a2
    ! Heun: a2 = 1/2
    ! Ralston: a2 = 2/3
    ! Half Point: a2 = 1
    ! ------------

    a2 = 1._dp/2._dp
    a1 = 1._dp - a2
    q = 1._dp/(2*a2)

    call update_ki(nequs, x(1), omega_0, phi0(1), omega_0, k1, func)
    call update_ki(nequs, x(1)+q*step, phi0(1)+q*step*k1(1), phi0, omega_0+q*step*k1(nequs), k2, func)

    do j = 1, nequs - 1
        phi(1, j) = phi0(j) + step*(a1*k1(j) + a2*k2(j))
    end do
    omega(1) = omega_0 + step*(a1*k1(nequs) + a2*k2(nequs))

    do i = 2, n_points
        call update_ki(nequs, x(i-1), omega(i-1), phi(i-1,:), omega(i-1), k1, func)
        call update_ki(nequs, x(i-1)+q*step, phi(i-1,1)+q*step*k1(1), phi(i-1,:), omega(i-1)+q*step*k1(nequs), k2, func)

        do j = 1, nequs - 1
            phi(i, j) = phi(i - 1, j) + step*(a1*k1(j) + a2*k2(j))
        end do
        omega(i) = omega(i - 1) + step*(a1*k1(nequs) + a2*k2(nequs))
    end do
end subroutine RungeKutta2order

subroutine RungeKutta4order(n_points, nequs, step, x, phi0, omega_0, phi, omega, func)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, nequs, i, j
    real(dp):: step, k1(nequs), k2(nequs), k3(nequs), k4(nequs), phi0(nequs - 1), omega_0, &
               phi(n_points, nequs - 1), omega(n_points), func, x(n_points) ! x for V = V(x)
    external :: func

    ! We got two eqs -> phi' = omega; alpha = 2·phi(x)·[V(x) - E] (see schrodinger_eq for deduction)
    
    ! --- RK4 ---
    ! k1 = f(x0, y0)
    ! k2 = f(x0 + h/2, y0 + h·k1/2)
    ! k3 = f(x0 + h/2, y0 + h·k2/2)
    ! k4 = f(x0 + h, y0 + h·k3)
    ! y1 = y0 + h(k1 + 2·k2 + 2·k3 + k4)/6
    ! ------------

    call update_ki(nequs, x(1), omega_0, phi0(1), omega_0, k1, func)
    call update_ki(nequs, x(1)+step/2._dp, phi0(1)+step*k1(1)/2._dp, phi0, omega_0+step*k1(nequs)/2._dp, k2, func)
    call update_ki(nequs, x(1)+step/2._dp, phi0(1)+step*k2(1)/2._dp, phi0, omega_0+step*k2(nequs)/2._dp, k3, func)
    call update_ki(nequs, x(1)+step, phi0(1)+step*k3(1), phi0, omega_0+step*k3(nequs), k4, func)

    do j = 1, nequs - 1
        phi(1, j) = phi0(j) + step*(k1(j) + 2._dp*k2(j) + 2._dp*k3(j) + k4(j))/6._dp
    end do
    omega(1) = omega_0 + step*(k1(nequs) + 2._dp*k2(nequs) + 2._dp*k3(nequs) + k4(nequs))/6._dp

    do i = 2, n_points
        call update_ki(nequs, x(i-1), omega(i-1), phi(i-1,:), omega(i-1), k1, func)
        call update_ki(nequs, x(i-1)+step/2._dp, phi(i-1,1)+step*k1(1)/2._dp, phi(i-1,:), omega(i-1)+step*k1(nequs)/2._dp, k2, func)
        call update_ki(nequs, x(i-1)+step/2._dp, phi(i-1,1)+step*k2(1)/2._dp, phi(i-1,:), omega(i-1)+step*k2(nequs)/2._dp, k3, func)
        call update_ki(nequs, x(i-1)+step, phi(i-1,1)+step*k3(1), phi(i-1,:), omega(i-1)+step*k3(nequs), k4, func)

        do j = 1, nequs - 1
            phi(i, j) = phi(i - 1, j) + step*(k1(j) + 2._dp*k2(j) + 2._dp*k3(j) + k4(j)) / 6._dp
        end do
        omega(i) = omega(i - 1) + step*(k1(nequs) + 2._dp*k2(nequs) + 2._dp*k3(nequs) + k4(nequs)) / 6._dp
    end do
end subroutine RungeKutta4order

subroutine Cash_Karp(n_points, nequs, step, x, phi0, omega_0, phi, omega, func) ! Wdoesn't work for some friking reason
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, nequs, i, j
    real(dp):: step, k1(nequs), k2(nequs), k3(nequs), k4(nequs), k5(nequs), k6(nequs), phi0(nequs - 1), omega_0, &
               phi(n_points, nequs - 1), omega(n_points), func, x(n_points), a(6), b(5, 6), c(6) ! x for V = V(x)
    external :: func, update_k

    ! We got two eqs -> phi' = omega; alpha = 2·phi(x)·[V(x) - E] (see schrodinger_eq for deduction)
    
    ! --- CK ---
    ! k1 = hf(x0, y0)
    ! ki = hf(x0 + a(i)·h, y0 + sum:J=1,i(b(j, i)·kj))
    ! y1 = y0 + h·sum:i=1,6(c(i)·ki)
    ! -----------

    a = [0._dp, 1._dp/5._dp, 3._dp/10._dp, 3._dp/5._dp, 1._dp, 7._dp/8._dp]
    b = reshape((/ 0._dp, 1._dp/5._dp, 3._dp/40._dp, 3._dp/10._dp, -11._dp/54._dp, 1631._dp/55296._dp, &
                   0._dp, 0._dp, 9._dp/40._dp, -9._dp/10._dp, 5._dp/2._dp, 175._dp/512._dp, &
                   0._dp, 0._dp, 0._dp, 6._dp/5._dp, -70._dp/27._dp, 575._dp/13824._dp, &
                   0._dp, 0._dp, 0._dp, 0._dp, 35._dp/27._dp, 44275._dp/110592._dp, &
                   0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 253._dp/4096._dp /), &
                   shape(b), order=(/ 2, 1 /))                   
    c = [37._dp/378._dp, 0._dp, 250._dp/621._dp, 125._dp/594._dp, 0._dp, 512._dp/1771._dp]

    call update_ki(nequs, x(1), omega_0, phi0, omega_0, k1, func)
    k1 = step*k1
    call update_ki(nequs, x(1)+a(2)*step, phi0(1)+b(1,2)*k1(1), phi0, omega_0+b(1,2)*k1(nequs), k2, func)
    k2 = step*k2
    call update_ki(nequs, x(1)+a(3)*step, phi0(1)+b(1,3)*k1(1)+b(2,3)*k2(1), phi0, &
                   omega_0+b(1,3)*k1(nequs)+b(2,3)*k2(nequs), k3, func)
    k3 = step*k3
    call update_ki(nequs, x(1)+a(4)*step, phi0(1)+b(1,4)*k1(1)+b(2,4)*k2(1)+b(3,4)*k3(1), phi0, &
                   omega_0+b(1,4)*k1(nequs)+b(2,4)*k2(nequs)+b(3,4)*k3(nequs), k4, func)
    k4 = step*k4
    call update_ki(nequs, x(1)+a(5)*step, phi0(1)+b(1,5)*k1(1)+b(2,5)*k2(1)+b(3,5)*k3(1)+b(4,5)*k4(1), phi0, &
                   omega_0+b(1,5)*k1(nequs)+b(2,5)*k2(nequs)+b(3,5)*k3(nequs)+b(4,5)*k4(1), k5, func)
    k5 = step*k5
    call update_ki(nequs, x(1)+a(6)*step, phi0(1)+b(1,6)*k1(1)+b(2,6)*k2(1)+b(3,6)*k3(1)+b(4,6)*k4(1)+b(5,6)*k5(1), phi0, &
                   omega_0+b(1,6)*k1(nequs)+b(2,6)*k2(nequs)+b(3,6)*k3(nequs)+b(4,6)*k4(1)+b(5,6)*k5(1), k6, func)
    k6 = step*k6

    do j = 1, nequs - 1
        phi(1, j) = phi0(j) + step*(c(1)*k1(j) + c(2)*k2(j) + c(3)*k3(j) + c(4)*k4(j) + c(5)*k5(j) + c(6)*k6(j))
    end do
    omega(1) = omega_0 + step*(c(1)*k1(nequs) + c(2)*k2(nequs) + c(3)*k3(nequs) + c(4)*k4(nequs) + c(5)*k5(nequs) + c(6)*k6(nequs))

    do i = 2, n_points
        call update_ki(nequs, x(1), omega(i-1), phi(i-1,:), omega(i-1), k1, func)
        k1 = step*k1
        call update_ki(nequs, x(1)+a(2)*step, phi(i-1,1)+b(1,2)*k1(1), phi(i-1,:), omega(i-1)+b(1,2)*k1(nequs), k2, func)
        k2 = step*k2
        call update_ki(nequs, x(1)+a(3)*step, phi(i-1,1)+b(1,3)*k1(1)+b(2,3)*k2(1), phi(i-1,:), &
                    omega(i-1)+b(1,3)*k1(nequs)+b(2,3)*k2(nequs), k3, func)
        k3 = step*k3
        call update_ki(nequs, x(1)+a(4)*step, phi(i-1,1)+b(1,4)*k1(1)+b(2,4)*k2(1)+b(3,4)*k3(1), phi(i-1,:), &
                    omega(i-1)+b(1,4)*k1(nequs)+b(2,4)*k2(nequs)+b(3,4)*k3(nequs), k4, func)
        k4 = step*k4
        call update_ki(nequs, x(1)+a(5)*step, phi(i-1,1)+b(1,5)*k1(1)+b(2,5)*k2(1)+b(3,5)*k3(1)+b(4,5)*k4(1), phi(i-1,:), &
                    omega(i-1)+b(1,5)*k1(nequs)+b(2,5)*k2(nequs)+b(3,5)*k3(nequs)+b(4,5)*k4(1), k5, func)
        k5 = step*k5
        call update_ki(nequs, x(1)+a(6)*step, phi(i-1,1)+b(1,6)*k1(1)+b(2,6)*k2(1)+b(3,6)*k3(1)+b(4,6)*k4(1)+b(5,6)*k5(1), &
                    phi(i-1,:), omega(i-1)+b(1,6)*k1(nequs)+b(2,6)*k2(nequs)+b(3,6)*k3(nequs)+b(4,6)*k4(1)+b(5,6)*k5(1), k6, func)
        k6 = step*k6

        do j = 1, nequs - 1
            phi(i, j) = phi(i - 1, j) + step*(c(1)*k1(j) + c(2)*k2(j) + c(3)*k3(j) + c(4)*k4(j) + c(5)*k5(j) + c(6)*k6(j))
        end do
        omega(i) = omega(i-1) + step*(c(1)*k1(nequs)+c(2)*k2(nequs)+c(3)*k3(nequs)+c(4)*k4(nequs)+c(5)*k5(nequs)+c(6)*k6(nequs))
    end do
end subroutine Cash_Karp

subroutine tir(n_points, nequs, x, E1, E2, phi0, omega_0, phi, omega, func, method)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/E, fact, beta

    integer :: nequs, n_points, n_iter, i
    real(dp) :: E1, E2, E3, phi0(nequs - 1), x(n_points), omega_0, phi(n_points, nequs - 1), omega(n_points), phiE1, phiE2, phiE3, &
                step, func, tol, E, fact, beta
    external :: method, func

    tol = 0.000001_dp

    n_iter = 100 ! safe guard just in case the method doesn't converge (should not happen)
    step = (x(n_points) - x(1))/ real(n_points - 1, dp)
    E = E1
    call method(n_points, nequs, step, x, phi0, omega_0, phi, omega, func)
    phiE1 = phi(n_points, 1)
    print *, phiE1
    
    E = E2
    call method(n_points, nequs, step, x, phi0, omega_0, phi, omega, func)
    phiE2 = phi(n_points, 1)

    E3 = (E1*phiE2 - E2*phiE1) / (phiE2 - phiE1)
    E = E3
    call method(n_points, nequs, step, x, phi0, omega_0, phi, omega, func)
    phiE3 = phi(n_points, 1)

    write(1, *) '# Eigenvalues'
    write(1, *) 1, E3, phiE3

    i = 1
    do while (abs(phiE3) .ge. tol)
        if (i .ge. n_iter) then
            exit
        end if

        phiE1 = phiE2
        phiE2 = phiE3
        E1 = E2
        E2 = E3
        E3 = (E1*phiE2 - E2*phiE1) / (phiE2 - phiE1)
        
        E = E3
        call method(n_points, nequs, step, x, phi0, omega_0, phi, omega, func)
        phiE3 = phi(n_points, 1)
        

        write(1, *) i+1, E3, phiE3
        i = i + 1
    end do
    E2 = E3
    write(1, *) ' '
    write(1, *) ' '
end subroutine tir

! ---------------------- EDO Presentation ----------------------
function schrodinger_eq(x, phi) result(alpha)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/E, fact, beta

    real(dp) :: phi, alpha, V, E, x, fact, beta

    ! Original: -phi''(x) / 2 + V·phi(x) = E·phi(x)
    ! Degree reduced: omega = phi'; omega' = phi'' -> -omega'(x) / 2 + V·phi(x) = E·phi(x) 
    ! omega'(t) = alpha = 2·phi(x)·[V(x) - E]

    alpha = (V(x) - E) * phi / fact ! p' = f(x, y)
end function schrodinger_eq

function schrodinger_eq_modified(x, phi) result(alpha)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/E, fact, beta

    real(dp) :: phi, V, E, x, fact, beta, alpha

    alpha = (V(x) + beta*sin(x) - E) * phi / fact ! p' = f(x, y)
end function schrodinger_eq_modified

function V(x) result(Vx)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    real(dp) :: x, Vx, alpha, delta, V0

    alpha = 2._dp
    delta = 0.4_dp
    V0 = -50._dp

    Vx = V0 * sinh(alpha) / (cosh(alpha) + cosh(x / delta))
end function V

subroutine simpson(n_points, step, y, AS)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, y(n_points), step
    integer :: n_points, i

    AS = y(1) + y(n_points)

    do i = 1, n_points - 1
        if (mod(i, 2) .eq. 0) then
            AS = AS + 2 * y(i)
        else
            AS = AS + 4 * y(i)
        end if
    end do
    AS = step * AS / 3._dp
end subroutine simpson