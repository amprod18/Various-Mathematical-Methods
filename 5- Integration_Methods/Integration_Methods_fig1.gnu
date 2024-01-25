set terminal png size 1920, 1080; set output 'Integration_Methods_fig1.png'
set autoscale


file = 'Integration_Methods_res.dat'

# --------------------- Absolute Error --------------------- 
set ylabel "Absolute Error of the Area"
set xlabel "Step Size"
set grid
set title "Integration Methods Error"
set key left top
set logscale xy

plot file index 0 using 1:2 with lp title "Trapeziod", \
     file index 0 using 1:3 with lp title "Simpson", \
     file index 0 using 1:4 with lp title "Simpson 3/8", \
     file index 0 using 1:5 with lp title "Boole",\
     file index 0 using 1:6 with lp title "Romberg"
