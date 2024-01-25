set terminal png size 640, 480; set output "EDO_Partials_D_fig4.png"

set title "Poisson Equation"
set xlabel "x (cm)"
set ylabel "y (cm)"
set grid
set autoscale xfix
set autoscale yfix
set autoscale cbfix


file = "EDO_Partials_D_res.dat"

plot file index 12 with image notitle

replot