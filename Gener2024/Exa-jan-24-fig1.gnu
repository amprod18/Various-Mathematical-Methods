set terminal png size 640, 480; set output "Exa-jan-24-fig1.png"

set title "Differential Equation Solution"
set xlabel "zeta(t)"
set ylabel "phi(t)"
set grid

file = "Exa-jan-24-res.dat"

plot file index 0 using 2:3 w lp title "Z0 = 0.3",\
     file index 2 using 2:3 w lp title "Z0 = 0.9"