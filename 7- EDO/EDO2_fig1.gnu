set terminal png size 640, 480; set output "EDO2_fig1.png"

set title "Schrödinger Equation Non-Solutions"


file = "EDO2_res.dat"

set ylabel "phi (A^{-1/2})"
set xlabel "x (A)"
set grid
set key bottom right

plot file index 0 using 1:2 title "PhiE1",\
     file index 1 using 1:2 title "PhiE2",\
     file index 2 using 1:2 title "PhiE3",\
     file index 3 using 1:2 title "PhiE4"