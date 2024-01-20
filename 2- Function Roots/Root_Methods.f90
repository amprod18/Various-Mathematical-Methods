program Root_Finding
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, n_points, n_iter_b, n_iter_r
    real(dp) :: a, epsilon, preci, Dmax_b, Dmax_r, xmax, ymax, E0, E1, Th, pi, E_D_max_b, E_D_max_r, &
                E_min, E_max, xmax_n, ymax_n, xmax_s, ymax_s
    real(dp), allocatable :: convergence_n(:), convergence_s(:), x(:), y(:), E(:), D(:), dD(:), t(:)
    external :: fun, fun2

    pi = acos(-1._dp)
    a = 189.857_dp
    epsilon = 0.995086_dp
    n_points = 120

    ! ------------------- Data Generation -------------------
    open(1, file="Root_Methods_res.dat")
    allocate(x(n_points), y(n_points), E(n_points), D(n_points), dD(n_points))

    write(1, *) '# Distance Function and Its derivate for Hale-Bopp Comet'
    do i = 1, n_points
        E(i) = (i - 1) * 2 * pi / (n_points - 1)
        call distance(a, epsilon, E(i), D(i), x(i), y(i))
    end do

    call numeric_derivate(n_points, E, D, dD)

    do i = 1, n_points
        write(1, *) E(i), D(i), dD(i)
    end do
    
    write(1, *) ' '
    write(1, *) ' '

    deallocate(x, y, E, D, dD)

    ! ------------------- Root Finding Methods -------------------
    
    preci = 1._dp/1000000000._dp
    print *, "Precision chosen: ", preci
    
    ! Bisection method used to find the maxima of the distance function
    E_min = 0.2_dp
    E_max = 6.1_dp
    n_iter_b = 1000
    n_iter_r = 1000

    write(1, *) '# Distance Maxima'
    call biseccio(fun, E_min, E_max, preci, n_iter_b, E_D_max_b)
    call distance(a, epsilon, E_D_max_b, Dmax_b, xmax, ymax)

    call regula_falsi(fun, E_min, E_max, preci, n_iter_r, E_D_max_r)
    call distance(a, epsilon, E_D_max_r, Dmax_r, xmax, ymax)

    print *, E_D_max_b, Dmax_b, n_iter_b
    print *, E_D_max_r, Dmax_r, n_iter_r

    write (1, *) E_D_max_b, Dmax_b, n_iter_b
    write (1, *) E_D_max_r, Dmax_r, n_iter_r

    write(1, *) ' '
    write(1, *) ' '

    ! Newton-Raphson method used to find the maxima of the distance function
    preci = 1._dp/100000000000._dp
    print *, "Precision chosen: ", preci
    n_points = 1000
    allocate(x(n_points), y(n_points), t(n_points), convergence_n(n_points), convergence_s(n_points))

    E0 = pi / 6._dp
    E1 = pi / 3._dp
    Th = 2526.5_dp

    write(1, *) '# Abnormal Excentricity and Orbit'
    do i = 1, n_points
        n_iter_r = 1000
        n_iter_b = 1000
        t(i) = (i - 1) * Th / (n_points - 1)

        ! Needs known derivate (in fun2)
        call newton_raphson(fun2, E0, preci, n_iter_r, convergence_n, Th, epsilon, t(i))
        call distance(a, epsilon, convergence_n(n_iter_r), Dmax_r, xmax_n, ymax_n)

        ! Doesn't need known derivate (in fun2)
        call secant_method(fun2, E0, E1, preci, n_iter_b, convergence_s, Th, epsilon, t(i))
        call distance(a, epsilon, convergence_s(n_iter_b), Dmax_b, xmax_s, ymax_s)

        write(1, *) t(i), convergence_n(n_iter_r), xmax_n, ymax_n, convergence_s(n_iter_b), xmax_s, ymax_s
    end do

    write(1, *) ' '
    write(1, *) ' '

    deallocate(convergence_n, convergence_s, x, y)
    close(1)

    call system('gnuplot -p Root_Methods_fig1.gnu')
    call system('gnuplot -p Root_Methods_fig2.gnu')
    call system('gnuplot -p Root_Methods_fig3.gnu')
    call system('gnuplot -p Root_Methods_fig4.gnu')
end program Root_Finding

subroutine fun(x, preci, valorf, valordf)
    ! Subroutine that returns the function and its corresponding derivate evaluated in the desired point
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp) :: x, preci, f, valorf, valordf

    valorf = f(x)
    valordf = (f(x + preci) - f(x - preci)) / (2._dp * preci)
end subroutine fun

function f(x) result(fres)
    ! Test function to apply the afforementioned root retrieval methods
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: fres, pi, epsilon

    pi = acos(-1._dp)
    epsilon = 0.995086
    fres = sin(2._dp * x)*(1._dp - epsilon**2) - (cos(x)*(2._dp - epsilon**2) - epsilon)*sin(x)
end function f

subroutine distance(a, epsilon, E, D, x, y)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), intent(in) :: a, epsilon, E
    real(dp), intent(out) ::  D
    real(dp) :: x, y

    x = a * (cos(E) - epsilon)
    y = a * sqrt(1._dp - epsilon**2)*sin(E)
    D = sqrt(x**2 + y**2)
end subroutine distance

function excentricity(t, Th, epsilon, E) result(fres)
    ! Test function to apply the afforementioned root retrieval methods
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), intent(in) :: Th, epsilon, t
    real(dp) :: pi, fres, E

    pi = acos(-1._dp)

    fres = E + epsilon*sin(E) - 2*pi*t / Th
end function excentricity

subroutine fun2(x, valorf, valordf, Th, epsilon, t)
    ! Subroutine that returns the function and its corresponding derivate evaluated in the desired point
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp) :: x, valorf, valordf, excentricity, Th, epsilon, t

    valorf = excentricity(t, Th, epsilon, x)
    valordf = 1 + epsilon * cos(x)
end subroutine fun2

subroutine numeric_derivate(n_points, x_values, funci, dfunci)
    ! Function that recieves 2 vectors of x values and the coresponding f(x) and returns the numeric derivate
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_points, i
    real(dp) :: x_values(n_points), funci(n_points), dfunci(n_points), h

    ! Assuming ordered x_values/funci pairs
    h = x_values(2) - x_values(1)
    dfunci(1) = (funci(2) - funci(1)) / h
    dfunci(n_points) = (funci(n_points) - funci(n_points - 1)) / h 

    do i = 2, n_points - 1
        dfunci(i) = (funci(i + 1) - funci(i - 1)) / (2._dp * h)
    end do
end subroutine numeric_derivate

subroutine biseccio(fun, A, B, preci, n_iter, root)
    ! Subroutine that implements the Bisection method to retrieve a solution of a funtion (f(x) = 0) 
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: A, B, preci, root, C, fa, fb, fc, diff, dfa, dfb, dfc
    integer :: n_iter, i
    external :: fun

    call fun(A, preci, fa, dfa)
    call fun(B, preci, fb, dfb)
    n_iter = ceiling(log((B - A) / preci) / log(2._dp))

    do i = 1, n_iter
        C = (A + B) / 2._dp
        call fun(C, preci, fc, dfc)
        diff = B - A
        if (fa * fb > 0) then
            print *, "Bisection method fails in the interval [", A, ", ", B, "]"
            return
        end if

        if (fc .EQ. 0.) then
            root = C
            return
        end if

        if (fa * fc < 0) then
            B = C
        else
            A = C
        end if

        if (diff < preci) then
            n_iter = i
            root = C
            exit
        end if
    end do
    root = C
end subroutine biseccio

subroutine regula_falsi(fun, A, B, preci, n_iter, root)
    ! Subroutine that implements the Bisection method to retrieve a solution of a funtion (f(x) = 0) 
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: A, B, preci, root, C, fa, fb, fc, diff, dfa, dfb, dfc
    integer :: n_iter, i
    external :: fun

    call fun(A, preci, fa, dfa)
    call fun(B, preci, fb, dfb)
    n_iter = ceiling(log((B - A) / preci) / log(2._dp))

    do i = 1, n_iter
        C = (A*fb + B*fa) / (fb - fa)
        call fun(C, preci, fc, dfc)
        diff = min(C - A, B - C)
        if (fa * fb > 0) then
            print *, "Regula-Falsi method fails in the interval [", A, ", ", B, "]"
            return
        end if

        if (fc .EQ. 0.) then
            root = C
            return
        end if

        if (fa * fc < 0) then
            B = C
        else
            A = C
        end if

        if (diff < preci) then
            n_iter = i
            root = C
            exit
        end if
    end do
    root = C
end subroutine regula_falsi

subroutine newton_raphson(fun, xini, preci, n_iter, root, Th, epsilon, t)
    ! Subroutine that implements the Newton-Raphson method to retrieve a solution of a funtion (f(x) = 0) 
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp) :: xini, preci, x, fx, dfx, error, Th, epsilon, t
    integer :: n_iter, i
    real(dp) :: root(n_iter)
    external :: fun

    x = xini
    do i = 1, n_iter
        call fun(x, fx, dfx, Th, epsilon, t)
        root(i) = x
        if (abs(fx / dfx) < preci) then
            exit
        end if
        x = x - fx / dfx
    end do
    n_iter = i
    error = abs(fx / dfx) ! Not in Use
end subroutine newton_raphson

subroutine secant_method(fun, xini1, xini2, preci, n_iter, root, Th, epsilon, t)
    ! Subroutine that implements the Newton-Raphson method to retrieve a solution of a funtion (f(x) = 0) 
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp) :: xini1, xini2, preci, xn, xn_1, fx, fx_1, dfx, error, Th, epsilon, t
    integer :: n_iter, i
    real(dp) :: root(n_iter)
    external :: fun

    xn = xini1
    xn_1 = xini2
    do i = 1, n_iter
        fx_1 = fx
        call fun(xn, fx, dfx, Th, epsilon, t)
        dfx = (xn - xn_1) / (fx - fx_1)
        root(i) = xn
        if (abs(fx / dfx) < preci) then
            exit
        end if
        xn_1 = xn
        xn = xn - fx / dfx
    end do
    n_iter = i
    error = abs(fx / dfx) ! Not in Use
end subroutine secant_method