set terminal png size 1920, 1080; set output "EDO_fig8.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set xlabel "t (s)"
set ylabel "phi (rad)"
set key bottom left
set grid

plot file index 1 using 1:2 title "Euler", \
     file index 1 using 1:4 title "Enhanced Euler",\
     file index 1 using 1:6 title "Adams-Bashforth",\
     file index 1 using 1:8 title "Adams 3p",\
     file index 1 using 1:10 title "Adams 4p",\
     file index 1 using 1:12 title "Adams-Moulton",\
     file index 1 using 1:14 title "Hamming"