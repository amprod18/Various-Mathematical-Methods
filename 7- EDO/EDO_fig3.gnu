set terminal png size 1920, 1080; set output "EDO_fig3.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set xlabel "phi (rad)"
set ylabel "omega (rad/s)"
set key top left
set grid

plot file index 1 using 2:3 title "Euler", \
     file index 1 using 4:5 title "Enhanced Euler",\
     file index 1 using 6:7 title "Adams-Bashforth",\
     file index 1 using 8:9 title "Adams 3p",\
     file index 1 using 10:11 title "Adams 4p",\
     file index 1 using 12:13 title "Adams-Moulton",\
     file index 1 using 14:15 title "Hamming"