set terminal png size 1920,1080; set output 'Interpolations_fig1.png'
set autoscale
set ylabel "Position (m)"
set xlabel "Time (s)"
set grid
set title "Pistons over time"

file = "Interpolation_res.dat"
plot for [i=2:6] file index 0 using 1:i w lp title "Piston ".(i-1)