set terminal png size 1920, 1080; set output "EDO_fig4.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set multiplot layout 2, 3

# plot 11
set ylabel "Energy (J) [phi0 = pi - 0.042]"
set title "Potential Energy (J)"
set grid

plot file index 2 using 1:2 title "Euler", \
     file index 2 using 1:3 title "Enhanced Euler",\
     file index 2 using 1:4 title "Adams-Bashforth"

# plot 12
set title "Kinetic Energy (J)"
set key top left
set grid

plot file index 2 using 1:5 title "Euler", \
     file index 2 using 1:6 title "Enhanced Euler", \
     file index 2 using 1:7 title "Adams-Bashforth"

# plot 13
set title "Total Energy (J)"
set key top left
set grid

plot file index 2 using 1:8 title "Euler", \
     file index 2 using 1:9 title "Enhanced Euler", \
     file index 2 using 1:10 title "Adams-Bashforth"

# plot 21
set xlabel "t (s)"
set ylabel "Energy (J) [phi0 = 1]"
set grid

plot file index 3 using 1:2 title "Euler", \
     file index 3 using 1:3 title "Enhanced Euler",\
     file index 3 using 1:4 title "Adams-Bashforth"

# plot 22
set xlabel "t (s)"
set key top left
set grid

plot file index 3 using 1:5 title "Euler", \
     file index 3 using 1:6 title "Enhanced Euler", \
     file index 3 using 1:7 title "Adams-Bashforth"

# plot 23
set xlabel "t (s)"
set key top left
set grid

plot file index 3 using 1:8 title "Euler", \
     file index 3 using 1:9 title "Enhanced Euler", \
     file index 3 using 1:10 title "Adams-Bashforth"

unset multiplot
replot