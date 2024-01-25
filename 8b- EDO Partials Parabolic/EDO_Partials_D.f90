program EDO_Partials_D
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/hx, ht, L, omega

    integer :: i, j, x_size, t_size, n_iter
    real(dp) :: L, hx, ht, T0(2), C, omega, heat_source_1, heat_source_2, initial_func
    real(dp), allocatable :: T(:, :)
    ! character(len=15) :: icontrol
    logical :: verbose, make_gif
    external :: heat_source_1, heat_source_2, initial_func

    ! Type of diff eqs
    ! dt^2 -> A; dxdt -> B; dx^2 -> C; dx,dt,x,t -> D 
    ! disc = B^2 - 4AC
    ! disc < 0 -> Eliptic , disc > 0 -> Hiperbolic, disc = 0 -> Parabolic

    n_iter = 1000
    L = 1._dp
    hx = 0.04_dp
    ht = 0.00075_dp
    C = 1._dp !kappa/c·rho 

    ! Dirichlet Conditions
    x_size = floor(L/hx) + 1 
    t_size = n_iter
    T0 = [0._dp, 0._dp] ! Borders: x=0, x=L; indepeendent of t

    open(1, file='EDO_Partials_D_res.dat')

    allocate(T(t_size, x_size))

    ! Poisson eq
    verbose = .False.
    make_gif = .False.

    write(1,*) '# Direct Method'
    call initiateT(T0, t_size, x_size, T, initial_func)
    call direct_method(t_size, x_size, T, C, verbose, make_gif, heat_source_1)
    do i = 1, t_size
        do j = 1, x_size
            write(1, *) (i-1)*ht, (j-1)*hx, T(i, j)
        end do
        write(1, *) ' '
    end do
    write(1,*) ' '

    write(1,*) '# Crank-Nicolson Method'
    call initiateT(T0, t_size, x_size, T, initial_func)
    call Crank_Nickolson(t_size, x_size, T, C, verbose, make_gif, heat_source_1)
    do i = 1, t_size
        do j = 1, x_size
            write(1, *) (i-1)*ht, (j-1)*hx, T(i, j)
        end do
        write(1, *) ' '
    end do
    write(1,*) ' '

    deallocate(T)
    
    close(1)
    call system('gnuplot -p EDO_Partials_D_fig1.gnu')
    call system('gnuplot -p EDO_Partials_D_fig2.gnu')
    ! call system('gnuplot -p EDO_Partials_D_fig3.gnu')
    ! call system('gnuplot -p EDO_Partials_D_fig4.gnu')
    ! call system('gnuplot -p EDO_Partials_D_fig5.gnu')
end program EDO_Partials_D

! ---------------------- EDO Solving Methods ----------------------
subroutine initiateT(T0, t_size, x_size, T, initial_func) ! Static Initial conditions (Dirichlet)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/h, Lx, Ly, omega

    integer :: i, t_size, x_size
    real(dp) :: h, Lx, Ly, T(t_size, x_size), T0(2), omega, initial_func

    T(1, 1) = T0(1)
    T(1, x_size) = T0(2)
    do i = 2, x_size - 1
        T(1, i) = initial_func(i) ! x-t0
    end do
end subroutine initiateT

subroutine direct_method(t_size, x_size, T, C, verbose, make_gif, rho)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/hx, ht, L, omega

    integer ::  i, j, k, max_iter, t_size, x_size
    real(dp) :: hx, ht, L, tol, max_var, T(t_size, x_size), newT(t_size, x_size), C, r, &
                omega, heat, rho
    character(len=15) :: iter
    logical :: verbose, make_gif
    external :: rho

    max_var = 0._dp
    tol = 0.001_dp
    max_iter = 1000
    newT = T
    r = C*ht/(hx**2)

    open(2, file='EDO_Partials_D_gif.dat')

    ! Iterate over interior, boundaries are fixed
    do i = 2, t_size
        do j = 2, x_size - 1
            heat = rho(j)
            newT(i, j) = r*(T(i-1, j+1) + T(i-1, j-1)) + (1-2*r)*T(i-1, j) + ht*heat

            if (abs(T(i, j) - newT(i, j)) .gt. max_var) then
                max_var = abs(T(i, j) - newT(i, j))
            end if
        end do
        T = newT
    end do
    
    ! Make a gif of the plate convergence
    if (make_gif) then
        write(iter, '(I5)') k
        write(2, *) '# Frame '//iter 

        do i = 1, t_size
            do j = 1, x_size
                write(2, *) (i-1)*ht, (j-1)*hx, T(i, j)
            end do
            write(2, *) ' '
        end do
        write(2, *) ' '
    end if

    ! Write a selected quadrant convergence
    if (verbose) then
        write(1, *) k, T(int(7.5_dp/ht), int(23.5_dp/hx))
    end if

    close(2)
    if (make_gif) then
        call system('gnuplot -p EDO_Partials_D_gif.gnu')
    end if
end subroutine direct_method

subroutine Crank_Nickolson(t_size, x_size, T, C, verbose, make_gif, rho)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/hx, ht, L, omega

    integer ::  i, j, k, max_iter, t_size, x_size
    real(dp) :: hx, ht, L, tol, max_var, T(t_size, x_size), newT(t_size, x_size), C, r, &
                omega, heat, rho
    character(len=15) :: iter
    logical :: verbose, make_gif
    external :: rho

    max_var = 0._dp
    tol = 0.001_dp
    max_iter = 1000
    newT = T
    r = C*ht/(hx**2)

    ! A*Ti+1,j = B
    ! A = (1 0 0 0...)(-r 2(1-r)-r -r 0 0 0...)(0 -r 2(1-r)-r -r 0 0 0...)...(0 0 0 ... 1) A()
    ! B = (T00 rTi,j-1+2(1-r)Ti,j+rTi,j+1 ... T0Lx) B(j)

    open(2, file='EDO_Partials_D_gif.dat')

    ! Iterate over interior, boundaries are fixed
    do i = 2, t_size
        do j = 2, x_size - 1
            heat = rho(j)
            newT(i, j) = r * T(i,j-1) + 2*(1._dp-r)*T(i,j) + r*T(i,j+1) + ht*heat

            if (abs(T(i, j) - newT(i, j)) .gt. max_var) then
                max_var = abs(T(i, j) - newT(i, j))
            end if
        end do
        T = newT
    end do
    
    ! Make a gif of the plate convergence
    if (make_gif) then
        write(iter, '(I5)') k
        write(2, *) '# Frame '//iter 

        do i = 1, t_size
            do j = 1, x_size
                write(2, *) (i-1)*ht, (j-1)*hx, T(i, j)
            end do
            write(2, *) ' '
        end do
        write(2, *) ' '
    end if

    ! Write a selected quadrant convergence
    if (verbose) then
        write(1, *) k, T(int(7.5_dp/ht), int(23.5_dp/hx))
    end if

    close(2)
    if (make_gif) then
        call system('gnuplot -p EDO_Partials_D_gif.gnu')
    end if
end subroutine Crank_Nickolson

! ---------------------- Heat Sources ----------------------
function heat_source_1(i) result(rhoxy)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/ht, hx, Lx, omega

    integer :: i
    real(dp) :: ht, hx, Lx, x, rhoxy, omega

    x = i*hx
    rhoxy = exp(-20*(x - 1._dp/2._dp)**2) - exp(-20*(x + 1._dp/2._dp)**2) - exp(-20*(x - 3._dp/2._dp)**2)
end function heat_source_1

function heat_source_2(i, j) result(rhoxy) ! No heatsource
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, j
    real(dp) :: rhoxy

    rhoxy = 0._dp
end function heat_source_2

function initial_func(index) result(T_out)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    common/dades/hx, ht, L, omega

    integer :: index
    real(dp) :: T_out, x, hx, ht, L, omega
    x = index*hx
    T_out = exp(-20*(x - 1._dp/2._dp)**2) - exp(-20*(x + 1._dp/2._dp)**2) - exp(-20*(x - 3._dp/2._dp)**2)
end function initial_func

