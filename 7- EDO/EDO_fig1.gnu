set terminal png size 1920, 1080; set output "EDO_fig1.png"

set title "Euler and Adams-Bashforth Methods"
set key top left

file = "EDO_res.dat"

set multiplot layout 2, 1

# Plot 11
set ylabel "phi (rad)"
set grid

plot file index 0 using 1:2 title "Euler", \
     file index 0 using 1:4 title "Enhanced Euler",\
     file index 0 using 1:6 title "Adams-Bashforth"

# Plot 21
set xlabel "t (s)"
set ylabel "omega (rad/s)"
set grid

plot file index 0 using 1:3 title "Euler", \
     file index 0 using 1:5 title "Enhanced Euler", \
     file index 0 using 1:5 title "Adams-Bashforth"

unset multiplot
replot