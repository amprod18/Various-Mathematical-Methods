set terminal gif animate delay 10 size 1920, 1080; set output "EDO_Partials_N.gif"

set title "Poisson Equation "
set xlabel "x (cm)"
set ylabel "y (cm)"
set grid
set autoscale xfix
set autoscale yfix
set autoscale cbfix

n_indices = 100
file = "EDO_Partials_N_gif.dat"

do for [i=0:n_indices] {plot file index i with image notitle}

replot