set terminal png size 1920, 1080; set output "figE3.png"

set title "RK4 Prevalence over Time"
set grid
set logscale y


file = "res.dat"

set multiplot layout 2, 3

# Plot 11
set ylabel "N = 10^4"
set title "p = 2/3"
plot file index 0 using 1:2 w lp title 'Simulation',\
     file index 1 using 1:2 w lp title 'RK4'

# Plot 12
set title "p = 1/2"
plot file index 2 using 1:2 w lp title 'Simulation',\
     file index 3 using 1:2 w lp title 'RK4'

# Plot 13
set title "p = 2/5"
plot file index 4 using 1:2 w lp title 'Simulation',\
     file index 5 using 1:2 w lp title 'RK4'

# Plot 21
set ylabel "N = 10^5"
set xlabel "t"
plot file index 6 using 1:2 w lp title 'Simulation',\
     file index 7 using 1:2 w lp title 'RK4'

# Plot 22
set xlabel "t"
plot file index 8 using 1:2 w lp title 'Simulation',\
     file index 9 using 1:2 w lp title 'RK4'

# Plot 23
set xlabel "t"
plot file index 10 using 1:2 w lp title 'Simulation',\
     file index 11 using 1:2 w lp title 'RK4'

unset multiplot
replot