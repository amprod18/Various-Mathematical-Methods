set terminal png size 640, 480; set output "EDO2_fig2.png"

set title "Schrödinger Equation Solutions"


file = "EDO2_res.dat"

set ylabel "E (eV)"
set xlabel "iteracion"
set grid

plot file index 4 using 1:2 with lp title "E1", \
     file index 6 using 1:2 with lp title "E2",\
     file index 8 using 1:2 with lp title "E3"