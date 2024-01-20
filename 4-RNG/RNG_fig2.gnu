set terminal png size 640, 480; set output "RNG_fig2.png"

set title "Montecarlo PDF vs Exact"
set xlabel "x (um^{-1})"
set ylabel "p(x)"
set grid
set yrange [0:0.12]

file = "RNG_res.dat"

plot file index 6 using 1:2 title "Exact distribution" with lp, \
     file index 5 using 1:2:4 with yerrorbars t "Montecarlo", \
     file index 5 using 1:2 with histeps t "PDF"

replot