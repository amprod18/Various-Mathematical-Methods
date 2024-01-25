set terminal png size 640, 480; set output "figE2.png"

set title "RK4 Prevalence over Time"
set xlabel "t"
set ylabel "Prevalence"
set grid
set logscale y


file = "res.dat"

plot file index 1 using 1:2 w lp title 'N=10^4;p=2/3',\
     file index 3 using 1:2 w lp title 'N=10^4;p=1/2',\
     file index 5 using 1:2 w lp title 'N=10^4;p=2/5',\
     file index 7 using 1:2 w lp title 'N=10^5;p=2/3',\
     file index 9 using 1:2 w lp title 'N=10^5;p=1/2',\
     file index 11 using 1:2 w lp title 'N=10^5;p=2/5'