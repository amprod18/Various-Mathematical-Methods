set terminal png size 640, 480; set output "EDO2_fig3.png"

set title "Normalized Schrödinger Equation Solutions"


file = "EDO2_res.dat"

set ylabel "phi (A^{-1/2})"
set xlabel "x (A)"
set grid
set key bottom right

plot file index 5 using 1:2 with lp title "Eigenvalue E1", \
     file index 7 using 1:2 with lp title "Eigenvalue E2",\
     file index 9 using 1:2 with lp title "Eigenvalue E3"