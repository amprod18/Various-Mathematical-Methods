set terminal png size 640, 480; set output 'Derivatives_fig3.png'
set autoscale

set ylabel "-exp(-x)"
set xlabel "x"
set grid
set title "Derivative with Different Methods"
set key left top

file = 'Derivatives_res.dat'

plot file index 2 using 1:2 with lp title "Analytic", \
     file index 2 using 1:3 with lp title "Numeric"