set terminal png size 1920,1080; set output 'Interpolations_fig2.png'
set autoscale
set ylabel "Position (m)"
set xlabel "Time (s)"
set grid
set title "Pistons over time"
set xrange [0:3]


file = "Interpolation_res.dat"

set multiplot layout 2, 3

do for [i=2:6] {
    plot file index 0 using 1:i w lp title "Piston ".(i-1),\
         file index i-1 using 1:2 w lp title "Linear Interpolation ".(i-1),\
         file index i-1 using 1:3 w lp title "Constant Interpolation ".(i-1),\
         file index i-1 using 1:4 w lp title "Polinomic Interpolation ".(i-1)
}