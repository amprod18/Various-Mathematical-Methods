set terminal png size 640, 480; set output "figE3.png"

set title "Inversion method vs Exact"
set xlabel "x"
set ylabel "p(x)"
set grid
set yrange [0:3.5]
set xrange [0:1.5]
set logscale y

file = "res1a.dat"

plot file index 3 using 5:6 w lp t "Exact PDF", \
     file index 3 using 1:2:4 with yerrorbars notitle, \
     file index 3 using 1:2 with boxes t "Montecarlo" lc rgb 'blue' smooth freq
replot