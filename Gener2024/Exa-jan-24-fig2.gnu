set terminal png size 640, 480; set output "Exa-jan-24-fig2.png"

set title "Montecarlo Integration and its Std"
set xlabel "N"
set ylabel "Integral / Error"
set grid

file = "Exa-jan-24-res.dat"

plot file index 4 using 1:2:3 w yerrorbars title "Integral MC",\
     file index 4 using 1:2 w lp notitle