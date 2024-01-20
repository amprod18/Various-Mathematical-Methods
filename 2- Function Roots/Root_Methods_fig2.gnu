set terminal png size 640, 480; set output 'Root_Methods_fig2.png'
set autoscale
set ylabel "y (UA)"
set xlabel "x (UA)"
set grid
set title "Orbit around the Sun of comet Hale-Bopp"

file = 'Root_Methods_res.dat'

plot file using 3:4 index 2 w lp title 'Newton-Raphson',\
     file using 6:7 index 2 w lp title 'Secant'