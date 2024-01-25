set terminal png size 640, 480; set output "EDO_Partials_N_fig2.png"

set autoscale xfix
set autoscale yfix
set autoscale cbfix
set key bottom left

file = "EDO_Partials_N_res.dat"

# First plot
set title 'T_0 = 120ºC'
set grid
set ylabel "T (ºC)"
set xlabel "N iter"
plot file index 4 with lp title 'Jacobi',\
     file index 5 with lp title 'Gauss-Seidel',\
     file index 6 with lp title 'Sobre-Relaxacio',\
     file index 7 with lp title '9P Laplacian'