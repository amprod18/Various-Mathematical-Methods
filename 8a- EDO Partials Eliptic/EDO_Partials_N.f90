program EDO_Partials_N
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/h, Lx, Ly, omega

    integer :: i, j, x_size, y_size
    real(dp) :: Lx, Ly, h, dT0(4), tinit(3), omega, heat_source_1, heat_source_2
    real(dp), allocatable :: T(:, :)
    character(len=15) :: icontrol
    logical :: verbose, make_gif
    external :: heat_source_1, heat_source_2

    ! Type of diff eqs
    ! dt^2 -> A; dxdt -> B; dx^2 -> C; dx,dt,x,t -> D 
    ! disc = B^2 - 4AC
    ! disc < 0 -> Eliptic , disc > 0 -> Hiperbolic, disc = 0 -> Parabolic

    Lx = 33.5_dp
    Ly = 45.5_dp
    h = 0.25_dp
    
    ! Neumann Conditions
    x_size = int(Lx/h) + 2
    y_size = int(Ly/h) + 2
    dT0 = [0.5_dp, 17._dp, -11.2_dp, -25.3_dp] ! Borders: x-0, 0-y, x-Ly, Lx-y

    tinit = [10._dp, 120._dp, 1040._dp]

    open(1, file='EDO_Partials_N_res.dat')

    allocate(T(x_size, y_size))

    ! Poisson eq
    verbose = .True.
    make_gif = .False.
    do i = 1, size(tinit)
        do j = 1, 4
            if (j .eq. 1) then
                icontrol = 'jacobi' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
                write(1,*) '# Jacobi Method'
                call initiateT(dT0, x_size, y_size, T, tinit(i))
                call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_1)
                write(1,*) ' '
                write(1,*) ' '
            else if (j .eq. 2) then
                icontrol = 'gauss_seidel' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
                write(1,*) '# Gauss-Seidel Method'
                call initiateT(dT0, x_size, y_size, T, tinit(i))
                call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_1)
                write(1,*) ' '
                write(1,*) ' '
            else if (j .eq. 3) then
                icontrol = 'sobre_relaxacio' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
                write(1,*) '# Sobre Relaxacio Method'
                omega = 1.35_dp
                call initiateT(dT0, x_size, y_size, T, tinit(i))
                call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_1)
                write(1,*) ' '
                write(1,*) ' '
            else
                icontrol = '9p_laplacian' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
                write(1,*) '# 9P Laplacian Method'
                omega = 1.35_dp
                call initiateT(dT0, x_size, y_size, T, tinit(i))
                call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_1)
                write(1,*) ' '
                write(1,*) ' '
            end if
        end do
    end do

    verbose = .False.
    make_gif = .False. ! make_gif = .True.
    write(1,*) '# Final Heat Map With Heatsources'
    icontrol = '9p_laplacian' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
    omega = 1.35_dp
    call initiateT(dT0, x_size, y_size, T, tinit(1))
    call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_1)
    do i = 2, x_size - 1
        do j = 2, y_size - 1
            write(1, *) (i-1)*h, (j-1)*h, T(i, j)
        end do
        write(1, *) ' '
    end do
    write(1,*) ' '

    make_gif = .False.
    write(1,*) '# Final Heat Map Without Heatsources'
    icontrol = '9p_laplacian' ! jacobi; gauss_seidel; sobre_relaxacio; 9p_laplacian
    call initiateT(dT0, x_size, y_size, T, tinit(1))
    call poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, heat_source_2)
    do i = 2, x_size - 1
        do j = 2, y_size - 1
            write(1, *) (i-1)*h, (j-1)*h, T(i, j)
        end do
        write(1, *) ' '
    end do
    write(1,*) ' '

    deallocate(T)
    
    close(1)
    call system('gnuplot -p EDO_Partials_N_fig1.gnu')
    call system('gnuplot -p EDO_Partials_N_fig2.gnu')
    call system('gnuplot -p EDO_Partials_N_fig3.gnu')
    call system('gnuplot -p EDO_Partials_N_fig4.gnu')
    call system('gnuplot -p EDO_Partials_N_fig5.gnu')
end program EDO_Partials_N

! ---------------------- EDO Solving Methods ----------------------
subroutine initiateT(T0, x_size, y_size, T, tinit) ! Static Initial conditions (Dirichlet)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/h, Lx, Ly, omega

    integer :: i, j, x_size, y_size
    real(dp) :: h, Lx, Ly, T(x_size, y_size), T0(4), tinit, omega

    do i = 2, x_size - 1
        do j = 2, y_size - 2
            T(i, j) = tinit
        end do
    end do

    do i = 1, x_size
        T(i, 1) = T0(1) ! x-0
    end do
    do j = 1, y_size
        T(1, j) = T0(2) ! 0-y
    end do
    do i = 1, x_size
        T(i, y_size) = T0(3) ! x-Ly
    end do
    do j = 1, y_size
        T(x_size, j) = T0(4) ! Lx-y
    end do
end subroutine initiateT

subroutine poisson_N(x_size, y_size, T, icontrol, verbose, make_gif, rho)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/h, Lx, Ly, omega

    integer ::  i, j, k, max_iter, x_size, y_size
    real(dp) :: h, Lx, Ly, tol, max_var, T(x_size, y_size), newT(x_size, y_size), &
                omega, heat, rho
    character(len=15) :: icontrol, iter
    logical :: verbose, make_gif
    external :: rho

    max_var = 0._dp
    tol = 0.1_dp
    max_iter = 1000
    newT = T

    open(2, file='EDO_Partials_N_gif.dat')

    do k = 1, max_iter
        ! Iterate over exterior boundaries
        do i = 2, x_size - 1
            newT(i, 2) = T(i, 3) + h * T(i, 1) ! x-0
        end do
        do j = 2, y_size - 1
            newT(2, j) = T(3, j) + h * T(1, j) ! 0-y
        end do
        do i = 2, x_size - 1
            newT(i, y_size) = T(i, y_size-2) + h * T(i, x_size) ! x-Ly
        end do
        do j = 2, y_size - 1
            newT(x_size, j) = T(x_size-2, j) + h * T(y_size, j) ! Lx-y
        end do
        ! Iterate over interior, boundaries have fixed flux (dT0)
        do i = 3, x_size - 2
            do j = 3, y_size - 2
                heat = rho(i, j)
                if (icontrol .eq. 'jacobi') then
                    newT(i, j) = (T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1) + h**2 * heat)/4._dp
                else if (icontrol .eq. 'gauss_seidel') then
                    newT(i, j) = (T(i+1, j) + newT(i-1, j) + T(i, j+1) + newT(i, j-1) + h**2 * heat)/4._dp
                else if (icontrol .eq. 'sobre_relaxacio') then
                    newT(i, j) = T(i, j) + &
                    omega * (T(i+1, j) + newT(i-1, j) + T(i, j+1) + newT(i, j-1) + h**2 * heat - 4._dp * T(i, j))/4._dp
                else if (icontrol .eq. '9p_laplacian') then
                    newT(i, j) = (T(i+1, j+1) + T(i-1, j-1) + T(i-1, j+1) + T(i+1, j-1) + &
                                  4*(T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1)) + 6 * h**2 * heat)/20._dp
                end if

                if (abs(T(i, j) - newT(i, j)) .gt. max_var) then
                    max_var = abs(T(i, j) - newT(i, j))
                end if
            end do
        end do
        T = newT
        
        ! Make a gif of the plate convergence
        if (make_gif) then
            write(iter, '(I5)') k
            write(2, *) '# Frame '//iter 

            do i = 1, x_size
                do j = 1, y_size
                    write(2, *) (i-1)*h, (j-1)*h, T(i, j)
                end do
                write(2, *) ' '
            end do
            write(2, *) ' '
        end if

        ! Write a selected quadrant convergence
        if (verbose) then
            write(1, *) k, T(int(7.5_dp/h), int(23.5_dp/h))
        end if
        
        if (max_var .lt. tol) exit
        max_var = 0._dp
    end do

    close(2)
    if (make_gif) then
        call system('gnuplot -p EDO_Partials_N_gif.gnu')
    end if
end subroutine poisson_N
! ---------------------- Heat Sources ----------------------
function heat_source_1(i, j) result(rhoxy)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/h, Lx, Ly, omega

    integer :: i, j
    real(dp) :: h, Lx, Ly, x, y, rhoxy, rho1, rho2, rho3, rho20, rho30, r, omega

    x = i * h
    y = j * h

    if (((x .ge. 18._dp) .and. (x .le. 22._dp)) .and. ((y .ge. 29._dp) .and. (y .le. 35._dp))) then
        rho1 = 3._dp
    else
        rho1 = 0._dp
    end if

    rho20 = 10._dp ! ºC/cm^2
    r = sqrt((x - 8._dp)**2 + (y - 22.5_dp)**2) ! distance to rho1 center
    rho2 = rho20 * exp(-((r - 5._dp)/0.3_dp)**2)

    rho30 = 6._dp ! ºC/cm^2
    r = sqrt((x - 22._dp)**2 + (y - 10.5_dp)**2) ! distance to rho1 center
    rho3 = rho30 * exp(-((r - 4._dp)/0.8_dp)**2)

    rhoxy = rho1 + rho2 + rho3
end function heat_source_1

function heat_source_2(i, j) result(rhoxy) ! No heatsource
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer :: i, j
    real(dp) :: rhoxy

    rhoxy = 0._dp
end function heat_source_2

