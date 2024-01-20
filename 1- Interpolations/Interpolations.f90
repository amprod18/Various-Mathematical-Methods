program Interpolations
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    
    integer :: i, j, k, n_points, n_interpo
    real(dp) :: omega_0, L
    real(dp), allocatable :: t(:), x(:, :), t_interpo(:), x_out_linear(:), x_out0(:), x_out_poli(:)
    omega_0 = 5._dp
    L = 18.5_dp
    k = 5
    n_points = 51
    n_interpo = 200

    ! ------------------- Data Generation -------------------
    allocate(t(n_points), x(n_points, k)) ! x = time x pistons
    open(1, file='Interpolation_res.dat')
    write(1, *) '# Analytic Positions'

    do i = 1, size(t)
        t(i) = real(i - 1, dp) * 5._dp / real(n_points - 1, dp)
        call posipisto(k, omega_0, L, t(i), x(i, :))
        write(1, *) t(i), x(i, :)
    end do

    write(1, *) ' '
    write(1, *) ' '
    
    ! ------------------- Interpolations -------------------
    allocate(t_interpo(n_interpo), x_out_linear(n_interpo), x_out0(n_interpo), x_out_poli(n_interpo))
    
    write(1, *) '# Interpolated Positions'
    n_points = 31
    do i = 1, k
        write(1, *) '# Piston kth'

        call linear_interpolation(n_points, n_interpo, t, x(:n_points, i), t_interpo, x_out_linear)
        call constant_interpolation(n_points, n_interpo, t, x(:n_points, i), t_interpo, x_out0)
        call polinomic_interpolation(n_points, n_interpo, t, x(:n_points, i), t_interpo, x_out_poli)
        do j = 1, n_interpo
            write(1, *) t_interpo(j), x_out_linear(j), x_out0(j), x_out_poli(j)
            x_out_linear(j) = 0
            x_out0(j) = 0
            x_out_poli(j) = 0
        end do
        write(1, *) ' '
        write(1, *) ' '
    end do
    
    deallocate(t, x, t_interpo, x_out_linear, x_out0, x_out_poli)
    close(1)
    
    call system('gnuplot -p Interpolations_fig1.gnu')
    call system('gnuplot -p Interpolations_fig2.gnu')
end program Interpolations

function radimano(L, k) result(radii_k)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    integer :: k
    real(dp) :: L, radii_k
    
    radii_k = (L / real(k, dp)) - 0.5_dp
end function radimano

function phi(k) result(phi_k)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    integer :: k
    real(dp) :: pi, phi_k

    pi = acos(-1._dp)
    phi_k = (real(k, dp) / 5._dp) ** 2 * pi
end function phi

subroutine posipisto(n_pistons, omega_0, L, t_in, x_out)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    integer :: k, n_pistons
    real(dp) :: omega_0, L, t_in, x_out(n_pistons), radii_k, radimano, phi_k, phi

    do k = 1, n_pistons
        phi_k = phi(k) ! Get piston phase
        radii_k = radimano(L, k) ! Get piston radii
        ! Get piston position
        x_out(k) = radii_k * cos(omega_0*t_in + phi_k) + sqrt(L**2 - (radii_k * sin(omega_0*t_in + phi_k))**2)
    end do
end subroutine posipisto

subroutine linear_interpolation(n_points, n_interpo, x_in, y_in, x_out, y_out)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_points, n_interpo, i, j
    real(dp) :: x_in(n_points), y_in(n_points), x_out(n_interpo), y_out(n_interpo), slope

    do i = 1, n_interpo
        x_out(i) = (x_in(n_points) - x_in(1)) * (i - 1)/real(n_interpo - 1, dp) + x_in(1)
    end do

    j = 1
    do i = 1, n_interpo
        do while (j  <= n_points)
            if (x_in(j) .ge. x_out(i)) then
                slope = (y_in(j) - y_in(j-1)) / (x_in(j) - x_in(j-1))
                y_out(i) = slope * (x_out(i) -  x_in(j-1)) + y_in(j-1)
                j = n_points + 1
            end if
            j = j + 1
        end do
        j = 1
    end do
end subroutine linear_interpolation

subroutine constant_interpolation(n_points, n_interpo, x_in, y_in, x_out, y_out)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_points, n_interpo, i, j
    real(dp) :: x_in(n_points), y_in(n_points), x_out(n_interpo), y_out(n_interpo)

    do i = 1, n_interpo
        x_out(i) = (x_in(n_points) - x_in(1)) * (i - 1)/real(n_interpo - 1, dp) + x_in(1)
    end do

    j = 1
    do i = 1, n_interpo
        do while (j <= n_points)
            if (x_in(j) .ge. x_out(i)) then
                y_out(i) = y_in(j - 1)
                j = n_points + 1
            end if
            j = j + 1
        end do
        j = 1
    end do
end subroutine constant_interpolation

subroutine polinomic_interpolation(n_points, n_interpo, x_in, y_in, x_out, y_out)
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer :: n_points, n_interpo, i, j, k, f
    real(dp) :: x_in(n_points), y_in(n_points), x_out(n_interpo), y_out(n_interpo), prod, step

    do i = 1, n_interpo
        x_out(i) = (x_in(n_points) - x_in(1)) * (i - 1)/real(n_interpo - 1, dp) + x_in(1)
    end do
    
    step = x_in(2) - x_in(1)
    do i = 1, n_interpo
        do while (j <= n_points)
            if (x_in(j) .ge. x_out(i)) then
                do f = j-1, j+1
                    prod = 1._dp
                    do k = j-1, j+1
                        if (k .ne. f) then
                            prod = prod * (x_out(i) - x_in(k)) / (x_in(f) - x_in(k))
                        end if
                    end do
                    y_out(i) = y_out(i) + prod * y_in(f)
                end do
                j = n_points + 1
            end if
            j = j + 1
        end do
        j = 1
    end do
end subroutine polinomic_interpolation

! Add more interpolations