program Integration_Methods
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AT, AS, AS38, AB, ARB, Areal, xini, xend, a, b, h, pi, YKohoutek
    integer :: i
    external :: YKohoutek

    a = 508.633_dp
    b = 429.074_dp
    pi = acos(-1._dp)
    xini = -4._dp * a
    xend = -7._dp * a/2._dp
    
    open(1, file="Integration_Methods_res.dat")

    Areal = a * b * (3._dp * sqrt(3._dp) + 2._dp * pi) / 24._dp
    print *, Areal
    do i = 3, 20, 2
        h = (xend - xini) / real(2 ** i, dp)
        call trapezi(xini, xend, i, AT, YKohoutek)
        call simpson(xini, xend, i, AS, YKohoutek)
        call simpsontresvuit(xini, xend, i, AS38, YKohoutek)
        call boole(xini, xend, i, AB, YKohoutek)
        call romberg(xini, xend, i, ARB, YKohoutek)

        write(1, *) h, abs(AT - Areal), abs(AS - Areal), abs(AS38 - Areal), abs(AB - Areal), &
                       abs(ARB - Areal), AT/Areal, AS/Areal, AS38/Areal, AB/Areal, ARB/Areal
                                   
    end do

    close(1)

    call system('gnuplot -p Integration_Methods_fig1.gnu')
    
end program Integration_Methods

! -------------------- Closed Integration Methods --------------------
subroutine trapezi(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, fa, fb, a, b, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp

    do i = 1, n_iter
        a = real(i - 1, dp) * step + xini
        b = real(i, dp) * step + xini
        ! print *, "Trapezi:", b

        fa = funcio(a)
        fb = funcio(b)

        AS = AS + step * (fa + fb) / 2._dp
    end do
end subroutine trapezi

subroutine simpson(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, a, fa, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = funcio(xini) + funcio(xend)

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

subroutine simpsontresvuit(xini, xend, m, AS, fun)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: i, m, k, k2
    real(dp) :: xini, xend, AS, fun, step
    
    k = 2**m
    k2 = ((k+2)/3)*3
    step = (xend-xini) / real(k2, dp)
    
    AS = fun(xini) + fun(xend)
    do i = 1, (k2 - 2), 3
        AS = AS + 3._dp*fun(xini + step*real(i, dp))
    end do
    do i = 2, (k2 - 1), 3
        AS = AS + 3._dp*fun(xini + step*real(i, dp))
    end do
    do i = 3, (k2 - 3), 3
        AS = AS + 2._dp*fun(xini + step*real(i, dp))
    end do
    AS = 3._dp*step*AS/8._dp
end subroutine simpsontresvuit

subroutine boole(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, a, fa, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = ((2**m + 3)/4)*4
    step = (xend - xini) / real(n_iter, dp)
    AS = 7._dp*(funcio(xini) + funcio(xend))

    do i = 1, n_iter - 1, 2
        a = real(i, dp) * step + xini
        fa = funcio(a)
        AS = AS + 32 * fa
    end do
    do i = 2, n_iter - 2, 4
        a = real(i, dp) * step + xini
        fa = funcio(a)
        AS = AS + 12 * fa
    end do
    do i = 4, n_iter - 4, 4
        a = real(i, dp) * step + xini
        fa = funcio(a)
        AS = AS + 14 * fa
    end do
    AS = 2 * step * AS / 45._dp
end subroutine boole

subroutine romberg(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, xini, xend, funcio, prev_col(m + 1), cur_col(m + 1)
    integer :: i, j, m
    external :: funcio
    
    do i = 0, m
        AS = 0._dp
        call trapezi(xini, xend, i+1, AS, funcio)
        prev_col(i+1) = AS
    end do

    do i = 1, m
        do j = 1, m - i + 1
            cur_col(j) = (4._dp**i * prev_col(j+1) - prev_col(j)) / real(4**i - 1, dp)
        end do

        do j = 1, m - i + 1
            prev_col(j) = cur_col(j)
        end do
    end do
    AS = prev_col(1)
end subroutine romberg

! -------------------- Open Integration Methods --------------------
subroutine gauss_legendre(xini, xend, m, AS, funcio)
    ! WTF
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, xini, xend, funcio, omega(2 * m - 1), zeta(2 * m - 1)
    integer :: i, m

    ! if (xini, xend) is not (-1, 1) make u = -1 + 2*((x-xini)/(xend-xini))
    omega = [1, 2, 3, 4] ! integral from -1 to 1 of Lagrange interpolation polinomia
    zeta = [1, 2, 3, 4] ! Zeros of the Legendre polinomia

    AS = 0._dp
    do i = 1, m
        AS = AS + omega(i) * funcio(zeta(i))
    end do
end subroutine gauss_legendre

! -------------------- Indefinite Integration Methods --------------------

! Convert infinite interval into finite:
! [0, +inf) -> t = 2 * atan(x) / pi t = [0, 1) 
! [0, +inf) -> t = (1 + x) / (1 - x) t = [-1, 1)

function YKohoutek(x) result(fres)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), intent(in) :: x
    real(dp) :: fres, a, b

    a = 508.633_dp
    b = 429.074_dp

    fres = b * sqrt(1._dp - ((x + 4._dp * a)**2) / a**2)
end function YKohoutek

function Lagrange_Polinomia(x, xi, m) result(prod)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: x, xi, prod, xk
    integer :: i, m

    prod = 1._dp
    do i = 1, m
        if (xi .ne. xk) then
            prod = prod * (x - xk) / (xi - xk)
        end if
    end do
end function Lagrange_Polinomia