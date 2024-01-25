set terminal png size 640, 480; set output "EDO2_fig4.png"

set title "Schrödinger Equation Custom Potential Solutions"


file = "EDO2_res.dat"

set ylabel "phi (A^{-1/2})"
set xlabel "x (A)"
set grid
set key top left


plot file index 10 using 1:2 with lp title "Beta = 0", \
     file index 12 using 1:2 with lp title "Beta = 5",\
     file index 14 using 1:2 with lp title "Beta = 15"