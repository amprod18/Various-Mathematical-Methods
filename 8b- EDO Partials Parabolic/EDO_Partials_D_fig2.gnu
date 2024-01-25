set terminal png size 640, 480; set output "EDO_Partials_D_fig2.png"

set title "Laplace Equation"
set xlabel "t (s)"
set ylabel "x (cm)"
set grid
set autoscale xfix
set autoscale yfix
set autoscale cbfix


file = "EDO_Partials_D_res.dat"

plot file index 1 with image notitle

replot