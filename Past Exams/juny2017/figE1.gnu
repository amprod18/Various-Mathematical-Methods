set terminal png size 640, 480; set output "figE1.png"

set title "Simulation Prevalence over Time"
set xlabel "t"
set ylabel "Prevalence"
set grid
set logscale y


file = "res.dat"

plot file index 0 using 1:2 w lp title 'N=10^4;p=2/3',\
     file index 2 using 1:2 w lp title 'N=10^4;p=1/2',\
     file index 4 using 1:2 w lp title 'N=10^4;p=2/5',\
     file index 6 using 1:2 w lp title 'N=10^5;p=2/3',\
     file index 8 using 1:2 w lp title 'N=10^5;p=1/2',\
     file index 10 using 1:2 w lp title 'N=10^5;p=2/5'