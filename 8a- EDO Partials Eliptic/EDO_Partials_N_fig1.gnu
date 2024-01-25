set terminal png size 640, 480; set output "EDO_Partials_N_fig1.png"

set autoscale xfix
set autoscale yfix
set autoscale cbfix
set key top left

file = "EDO_Partials_N_res.dat"

# First plot
set title 'T_0 = 10ºC'
set grid
set ylabel "T (ºC)"
set xlabel "N iter"
plot file index 0 with lp title 'Jacobi',\
     file index 1 with lp title 'Gauss-Seidel',\
     file index 2 with lp title 'Sobre-Relaxacio', \
     file index 3 with lp title '9P Laplacian'