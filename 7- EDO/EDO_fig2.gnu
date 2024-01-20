set terminal png size 1920, 1080; set output "EDO_fig2.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set multiplot layout 2, 1

set xlabel "t (s)"
set ylabel "phi (rad)"
set key bottom left
set grid

plot file index 1 using 1:2 title "Euler", \
     file index 1 using 1:4 title "Enhanced Euler",\
     file index 1 using 1:6 title "Adams-Bashforth"

set xlabel "t (s)"
set ylabel "omega (rad/s)"
set key top right
set grid

plot file index 1 using 1:3 title "Euler", \
     file index 1 using 1:5 title "Enhanced Euler", \
     file index 1 using 1:7 title "Adams-Bashforth"

unset multiplot

replot