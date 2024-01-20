set terminal png size 1920, 1080; set output "EDO_fig8b.png"

set title "Euler and Adams-Bashforth Methods"


file = "EDO_res.dat"

set xlabel "t (s)"
set ylabel "omega (rad/s)"
set key top right
set grid

plot file index 1 using 1:3 title "Euler", \
     file index 1 using 1:5 title "Enhanced Euler", \
     file index 1 using 1:7 title "Adams-Bashforth",\
     file index 1 using 1:9 title "Adams 3p",\
     file index 1 using 1:11 title "Adams 4p",\
     file index 1 using 1:13 title "Adams-Moulton",\
     file index 1 using 1:15 title "Hamming"