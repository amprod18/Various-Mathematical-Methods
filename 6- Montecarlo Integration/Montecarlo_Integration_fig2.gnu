set terminal png size 640, 480; set output "Montecarlo_Integration_fig2.png"

set title "Montecarlo Integration Method Convergence"
set xlabel "N"
set ylabel "Stardard Deviation"
set grid

file = "Montecarlo_Integration_res.dat"

plot file index 1 using 1:3 w lp title "Std Integral 2"
     # file index 1 using 1:4 w lp title "Error Integral 2"