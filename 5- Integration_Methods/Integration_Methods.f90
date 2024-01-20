program Integration_Methods
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AT, AS, AS38, AB, ART, ARS, ARB, Areal, xini, xend, a, b, h, pi, YKohoutek
    integer :: i
    external :: YKohoutek

    a = 508.633_dp
    b = 429.074_dp
    pi = acos(-1._dp)
    xini = -4._dp * a
    xend = -7._dp * a/2._dp
    
    open(1, file="Integration_Methods_res.dat")

    Areal = a * b * (3._dp * sqrt(3._dp) + 2._dp * pi) / 24._dp
    do i = 2, 20, 2
        h = (xend - xini) / real(2 ** i, dp)
        call trapezi(xini, xend, i, AT, YKohoutek)
        call simpson(xini, xend, i, AS, YKohoutek)
        call simpsontresvuit(xini, xend, i, AS38, YKohoutek)
        call boole(xini, xend, i, AB, YKohoutek)
        call repeated_trapezoids(xini, xend, i, ART, YKohoutek)
        call repeated_simpson(xini, xend, i, ARS, YKohoutek)
        call romberg(xini, xend, i, ARB, YKohoutek)

        write(1, "(15f20.12)") h, abs(AT - Areal), abs(AS - Areal), abs(AS38 - Areal), abs(AB - Areal), abs(ART - Areal), &
                                  abs(ARS - Areal), abs(ARB - Areal), AT/Areal, AS/Areal, AS38/Areal, AB/Areal, ART/Areal, &
                                  ARS/Areal, ARB/Areal 
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
        b = a + step
        ! print *, "Trapezi:", b

        fa = funcio(a)
        fb = funcio(b)

        AS = AS + step * (fa + fb) / 2._dp
    end do
end subroutine trapezi

subroutine simpson(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, fa, fb, fc, a, b, c, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp

    do i = 1, n_iter - 1, 2
        a = real(i - 1, dp) * step + xini
        b = a + step
        c = a + 2 * step
        ! print *, "Simpson:", c

        fa = funcio(a)
        fb = funcio(b)
        fc = funcio(c)

        AS = AS + step * (fa + 4*fb + fc) / 3._dp
    end do
end subroutine simpson

subroutine simpsontresvuit(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, fa, fb, fc, fd, a, b, c, d, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp

    do i = 1, n_iter - 2, 3
        a = real(i - 1, dp) * step + xini
        b = a + step
        c = a + 2 * step
        d = c + 3 * step
        ! print *, "3/8:", d

        fa = funcio(a)
        fb = funcio(b)
        fc = funcio(c)
        fd = funcio(d)

        AS = AS + 3 * step * (fa + 3*(fb + fc) + fd) / 8._dp
    end do
end subroutine simpsontresvuit

subroutine boole(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, fa, fb, fc, fd, fe, a, b, c, d, e, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp

    do i = 1, n_iter - 4, 4
        a = real(i - 1, dp) * step + xini
        b = a + step
        c = a + 2 * step
        d = c + 3 * step
        e = c + 4 * step
        ! print *, "Boole:", e

        fa = funcio(a)
        fb = funcio(b)
        fc = funcio(c)
        fd = funcio(d)
        fe = funcio(e)

        AS = AS + 2 * step * (7*(fa + fe) + 32*(fb + fd) + 12*fc) / 45._dp
    end do
end subroutine boole

! -------------------- Repeated Integration Methods --------------------
subroutine repeated_trapezoids(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp
    AS = AS + (funcio(xini) + funcio(xend))/2._dp

    do i = 2, n_iter - 2
        AS = AS + funcio(xini + real(i - 1, dp) * step)
    end do
    AS = step * AS
end subroutine repeated_trapezoids

subroutine repeated_simpson(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, step, xini, xend, funcio
    integer :: n_iter, i, m

    n_iter = 2 ** m + 1 ! Needs to be even
    step = (xend - xini) / real(n_iter, dp)
    AS = 0._dp
    AS = AS + funcio(xini) + funcio(xend)

    do i = 2, n_iter - 2
        if (mod(i, 2) .eq. 0) then
            AS = AS + 4._dp * funcio(xini + real(i - 1, dp) * step)
        else
            AS = AS + 2._dp * funcio(xini + real(i - 1, dp) * step)
        end if
    end do
    AS = step * AS / 3._dp
end subroutine repeated_simpson

subroutine romberg(xini, xend, m, AS, funcio)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: AS, xini, xend, funcio, prev_col(m + 1), cur_col(m + 1)
    integer :: i, j, m
    external :: funcio

    AS = 0._dp

    do i = 0, m
        call trapezi(xini, xend, i, AS, funcio)
        prev_col(i+1) = AS
    end do

    do i = 1, m
        do j = 1, m - i + 1
            cur_col(j) = (4**i * prev_col(j+1) - prev_col(j)) / (4**i - 1)
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