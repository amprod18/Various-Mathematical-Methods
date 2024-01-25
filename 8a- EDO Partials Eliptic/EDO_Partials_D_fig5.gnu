set terminal png size 640, 480; set output "EDO_Partials_D_fig5.png"

set title "Laplace Equation"
set xlabel "x (cm)"
set ylabel "y (cm)"
set grid
set autoscale xfix
set autoscale yfix
set autoscale cbfix


file = "EDO_Partials_D_res.dat"

plot file index 13 with image notitle

replot