set terminal png size 1920, 1080; set output "EDO_fig6.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set multiplot layout 2, 2

# Plot 11
set title '300 Steps'
set ylabel "Total Energy (J)"
set grid

plot file index 5 using 1:4 notitle

# Plot 12
set title '1000 Steps'
set grid

plot file index 6 using 1:4 notitle

# Plot 21
set title '2200 Steps'
set xlabel "t (s)"
set ylabel "Total Energy (J)"
set grid

plot file index 7 using 1:4 notitle

# Plot 22
set title '14500 Steps'
set xlabel "t (s)"
set grid

plot file index 8 using 1:4 notitle

unset multiplot
replot