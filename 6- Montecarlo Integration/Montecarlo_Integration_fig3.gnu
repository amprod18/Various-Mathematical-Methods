set terminal png size 640, 480; set output "Montecarlo_Integration_fig3.png"

set title "Montecarlo Integration Method Convergence"
set xlabel "N"
set ylabel "Stardard Deviation"
set grid
set logscale xy

file = "Montecarlo_Integration_res.dat"

plot file index 2 using 1:2:($2-$3):($2+$3) w yerrorbars title "Integral 3 AR",\
     file index 2 using 1:2 w lines notitle,\
     file index 2 using 1:4:($4-$5):($4+$5) w yerrorbars title "Integral 3 M",\
     file index 2 using 1:4 w lines notitle