program RNG
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer :: i, iseed, ncaixes, n_points, m
    real(dp) :: px, x, xlow, xhigh, L, cotasup, funcio_electro, mean_e, median_e, std_dev_e, variance_e, boxsize, &
                pi, AS, AS_t, sigma, funcio_rubidi, mean_rb, median_rb, std_dev_rb, variance_rb
    real(dp) :: numeros_e(40000), numeros_rb(10000)
    real(dp), allocatable :: xhis(:), vhis(:), errhis(:)
    external :: funcio_electro, funcio_rubidi
                
    ! ------------ Electron Atom ------------
    ncaixes = 80
    cotasup = 1._dp
    L = 3._dp ! L is in nanometers so 3e-9
    n_points = 40000
    xlow = -3._dp * L
    xhigh = 3._dp * L
    pi = acos(-1._dp)

    iseed = 20347972
    call srand(iseed)

    ! --- Data Generation ---
    open(1, file='RNG_res.dat')
    write (1, *) '#Ex 1'

    call accept_reject(n_points, numeros_e, xlow, xhigh, cotasup, funcio_electro)

    ! --- Electron Statistics ---
    call mean(n_points, numeros_e, mean_e)
    call median(n_points, numeros_e, median_e)
    call std_dev(n_points, numeros_e, std_dev_e)
    call variance(n_points, numeros_e, variance_e)

    write (1, *) '# Mean - Median - Standard Deviation - Variance of the Electron'
    write (1, *) mean_e, median_e, std_dev_e, variance_e
    write (1, *) ' '
    write (1, *) ' '

    print *, 'Montecarlo Values ar method:  ', 'Mean:', mean_e, 'Median:', median_e, 'Standard Deviation:', std_dev_e, &
             'Variance:', variance_e

    ! --- Histogram Generation ---
    write(1, *) '# Histogram'
    allocate(xhis(ncaixes), vhis(ncaixes), errhis(ncaixes))
    call histograma(n_points, numeros_e, xlow, xhigh, ncaixes, xhis, vhis, errhis, boxsize)
    
    do i = 1, ncaixes
        write (1, *) xhis(i), vhis(i), boxsize, errhis(i)
    end do
    write (1, *) ' '
    write (1, *) ' '
    deallocate(xhis, vhis, errhis)

    ! --- Electron PDF ---
    write (1, *) '# Exact PDF'
    do i = 1, n_points
        x = (i - 1) * 6._dp * L / (n_points - 1) - 3._dp * L
        px = funcio_electro(x)
        write (1, *) x, px
    end do
    write (1, *) ' '
    write (1, *) ' '


    m = 10
    AS_t = 0._dp
    call simpson(xlow, xhigh, m, AS_t, funcio_electro)
    print *, 'Total area under electron PDF (simpson):', AS_t

    AS = 0._dp
    xhigh = L / 2._dp
    call simpson(xlow, xhigh, m, AS, funcio_electro)
    print *, 'P(-3L =< x =< L/2) =', AS 
    write (1, *) '# Normalization constant', 'P(-3L =< x =< L/2)'
    write (1, *) AS_t, AS
    write (1, *) ' '
    write (1, *) ' '

    ! ------------ Rubidium Atom ------------
    ncaixes = 70
    cotasup = 1._dp
    sigma = 4._dp ! sigma is in micrometers so 4e-6
    n_points = 10000
    xlow = -3._dp * sigma
    xhigh = 3._dp * sigma
    
    ! --- Data Generation ---
    call accept_reject(n_points, numeros_rb, xlow, xhigh, cotasup, funcio_rubidi)
    
    ! --- Electron Statistics ---
    call mean(n_points, numeros_rb, mean_rb)
    call median(n_points, numeros_rb, median_rb)
    call std_dev(n_points, numeros_rb, std_dev_rb)
    call variance(n_points, numeros_rb, variance_rb)
    
    write (1, *) '#Ex 2'
    write (1, *) '# Mean - Median - Standard Deviation - Variance of the Electron'
    write (1, *) mean_rb, median_rb, std_dev_rb, variance_rb
    write (1, *) ' '
    write (1, *) ' '

    print *, 'Exact Values ar method:', 'Mean:', 0, 'Variance:', sigma ** 2, 'Standard Deviation:', sigma
    print *, 'Montecarlo Values ar method:  ', 'Mean:', mean_rb, 'Median:', median_rb, 'Standard Deviation:', std_dev_rb, &
             'Variance:', variance_rb

    ! --- Histogram Generation ---
    write(1, *) '# Histogram'
    allocate(xhis(ncaixes), vhis(ncaixes), errhis(ncaixes))
    call histograma(n_points, numeros_rb, xlow, xhigh, ncaixes, xhis, vhis, errhis, boxsize)
    do i = 1, ncaixes
        write (1, *) xhis(i), vhis(i), boxsize, errhis(i)
    end do
    write (1, *) ' '
    write (1, *) ' '
    deallocate(xhis, vhis, errhis)

    ! --- Rubidium PDF ---
    write (1, *) '# Exact PDF'
    do i = 1, n_points
        x = (i - 1) * 6._dp * sigma / (n_points - 1) - 3._dp * sigma
        px = funcio_rubidi(x)
        write (1, *) x, px
    end do
    write (1, *) ' '
    write (1, *) ' '

    m = 10
    AS_t = 0._dp
    call simpson(xlow, xhigh, m, AS_t, funcio_rubidi)
    print *, 'Total area under rubidium PDF (simpson):', AS_t

    close(1)

    call system('gnuplot -p RNG_fig1.gnu')
    call system('gnuplot -p RNG_fig2.gnu')
end program RNG

subroutine histograma(ndat, xdat, xa, xb, nbox, xhis, vhis, errhis, boxsize)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: ndat, nbox, ierr, i, ibox, icount
    real(dp) :: xa, xb, boxsize
    real(dp) :: xdat(ndat), xhis(nbox), vhis(nbox), errhis(nbox)

    if ( xa .ge. xb ) then
        ierr = 1
        return
    end if

    boxsize = (xb - xa) / nbox

    icount = 0

    do i = 1, nbox
        vhis(i) = 0
        errhis(i) = 0
    end do

    do i = 1, ndat
        ! data between xa and xb
        if ( xdat(i) .ge. xa .and. xdat(i).le.xb ) then
            ibox = int((xdat(i) - xa) / boxsize) + 1 
            if (ibox .eq. nbox + 1) ibox = nbox ! puts xb into the last box
            vhis(ibox) = vhis(ibox) + 1
            icount = icount + 1
        end if
    end do
    
    if ( icount .eq. 0 ) then
        ierr = 2
        return
    end if

    ierr = 0
    print *, "Accepted: ", icount, "out of:", ndat

    do i = 1, nbox
        xhis(i) = xa + boxsize/2._dp + (i - 1)*boxsize ! central value
        errhis(i) = sqrt(vhis(i)/icount * (1._dp - vhis(i)/icount))/boxsize / sqrt(icount * 1._dp) ! errorbar
        vhis(i) = vhis(i) / icount / boxsize ! normalized value
    end do
end subroutine histograma

! -------------- Statistical Toolkit --------------
subroutine mean(n_points, x_values, mean_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points
    real(dp) :: x_values(n_points), mean_val

    mean_val = sum(x_values) / real(n_points, dp)
 end subroutine mean

 subroutine median(n_points, x_values, median_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i, j
    real(dp) :: temp_array(n_points), x_values(n_points), median_val

    temp_array = x_values

    do i = 1, n_points
       do j = 1, n_points-i
          if (temp_array(j) > temp_array(j+1)) then
             temp_array(j) = temp_array(j+1)
             temp_array(j+1) = temp_array(j)
          end if
        end do
    end do

    if (mod(n_points, 2) == 0) then
        median_val = (temp_array(n_points/2) + temp_array(n_points/2 + 1)) / 2.0_dp
    else
        median_val = temp_array(n_points/2 + 1)
    end if
end subroutine median

subroutine variance(n_points, x_values, variance_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: n_points
    real(dp) :: x_values(n_points), mean_val, variance_val
    
    call mean(n_points, x_values, mean_val)
    
    variance_val = sum((x_values - mean_val)**2) / (real(n_points, dp) - 1.0_dp)    
end subroutine variance

subroutine std_dev(n_points, x_values, std_dev_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    
    integer :: n_points
    real(dp) :: x_values(n_points), std_dev_val, variance_val
    
    call variance(n_points, x_values, variance_val)
    
    std_dev_val = sqrt(variance_val)
end subroutine std_dev

subroutine covariance_matrix(n_points, n_vars, x_values, covariance_matrix_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, n_vars, i, j, k
    real(dp) :: x_values(n_points, n_vars), mean_vali, mean_valj, covariance_matrix_val(n_vars, n_vars), sum_x

    ! Compute the covariance matrix
    covariance_matrix_val = 0.0_dp
    do i = 1, n_vars
        call mean(n_points, x_values(:, i), mean_vali)
        do j = i, n_vars
            call mean(n_points, x_values(:, j), mean_valj)
            sum_x = 0.0_dp
            do k = 1, n_points
                sum_x = sum_x + (x_values(k, i) - mean_vali) * (x_values(k, j) - mean_valj)
            end do
            covariance_matrix_val(i, j) = sum_x / real(n_points - 1, dp)
            covariance_matrix_val(j, i) = covariance_matrix_val(i, j)
        end do
    end do
end subroutine covariance_matrix

subroutine correlation_matrix(n_points, n_vars, x_values, correlation_matrix_val)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, n_vars, i, j
    real(dp) :: x_values(n_points, n_vars), std_dev_vali, std_dev_valj, &
                correlation_matrix_val(n_vars, n_vars), covariance_matrix_val(n_vars, n_vars)

    ! Compute the covariance matrix
    call covariance_matrix(n_points, n_vars, x_values, covariance_matrix_val)
    correlation_matrix_val = 0.0_dp
    do i = 1, n_vars
        call std_dev(n_points, x_values(:, i), std_dev_vali)
        do j = i, n_vars
            call std_dev(n_points, x_values(:, j), std_dev_valj)
            correlation_matrix_val(i, j) = covariance_matrix_val(i, j) / (std_dev_vali * std_dev_valj)
            correlation_matrix_val(j, i) = correlation_matrix_val(i, j)
        end do
    end do
end subroutine correlation_matrix

! -------------- Generation Methods --------------
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

subroutine inversion(n_points, x_out, xlow, xhigh, cotasup, funcio) ! TBD
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: n_points, i, j
    real(dp) :: xlow, xhigh, cotasup, funcio, px, x_out(n_points), x1, x2, p, x
    external :: funcio
end subroutine inversion

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

function funcio_electro(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    real(dp) :: px, x, L
    L = 3._dp

    px = ((x / L)**2) * (9._dp - (x/L)**2) * 5._dp/(324._dp * L)
end function funcio_electro

function funcio_rubidi(x) result(px)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    real(dp) :: px, x, sigma, pi
    pi = acos(-1._dp)
    sigma = 4._dp

    px = exp(-(x**2) / (2._dp*(sigma**2)))/sqrt(2._dp*pi*(sigma**2))
end function funcio_rubidi

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