set terminal png size 1920, 1080; set output "EDO_fig5.png"

set title "Transition zone with Adams-Bashforth Method"


file = "EDO_res.dat"

set xlabel "phi (rad)"
set ylabel "omega (rad/s)"
set key bottom right
set grid

plot file index 4 using 2:3 title "omega_{0}^{-}", \
     file index 4 using 4:5 title "omega_{0}^{+}",\