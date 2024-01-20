set terminal png size 640, 480; set output 'Derivatives_fig2.png'
set autoscale

set ylabel "-exp(-x)"
set xlabel "x"
set grid
set title "Derivative with Different Methods"
set key left top

file = 'Derivatives_res.dat'

plot file index 1 using 1:2 with lp title "Analytic", \
     file index 1 using 1:3 with lp title "Forward",\
     file index 1 using 1:4 with lp title "Backward",\
     file index 1 using 1:5 with lp title "Central",\
     file index 1 using 1:6 with lp title "5 Points"