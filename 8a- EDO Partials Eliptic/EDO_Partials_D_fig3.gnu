set terminal png size 640, 480; set output "EDO_Partials_D_fig3.png"

set autoscale xfix
set autoscale yfix
set autoscale cbfix
set key bottom left

file = "EDO_Partials_D_res.dat"

# First plot
set title 'T_0 = 1040ºC'
set grid
set ylabel "T (ºC)"
set xlabel "N iter"
plot file index 8 with lp title 'Jacobi',\
     file index 9 with lp title 'Gauss-Seidel',\
     file index 10 with lp title 'Sobre-Relaxacio',\
     file index 11 with lp title '9P Laplacian'