set terminal png size 640, 480; set output "RNG_fig1.png"

set title "Montecarlo PDF vs Exact"
set xlabel "x (nm^{-1})"
set ylabel "p(x)"
set grid
set yrange [0:0.12]

file = "RNG_res.dat"

plot file index 2 using 1:2 title "Exact distribution" with lp, \
     file index 1 using 1:2:4 with yerrorbars notitle, \
     file index 1 using 1:2 with boxes t "Montecarlo" lc rgb 'blue' smooth freq

replot