set terminal png size 640, 480; set output 'Derivatives_fig1.png'
set autoscale

set ylabel "exp(-x) / -exp(-x)"
set xlabel "x"
set grid
set title "Analytic function and its Derivative"
set key right top

file = 'Derivatives_res.dat'

plot file index 0 using 1:2 with lp title "f(x)=x·exp(-x)", \
     file index 0 using 1:3 with lp title "f'(x)=exp(-x)-x·exp(-x)",\
     file index 0 using 1:4 with lp title "f''(x)=-2exp(-x)+x·exp(-x)"