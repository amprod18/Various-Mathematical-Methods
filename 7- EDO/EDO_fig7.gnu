set terminal png size 1920, 1080; set output "EDO_fig7.png"

set title "System Energy dependant on Step Size"

file = "EDO_res.dat"

set ylabel "Total Energy (J)"
set xlabel 't (s)'
set key top left
set grid

plot file index 5 using 1:3 w lp title '300 Steps', \
     file index 6 using 1:3 w lp title '1000 Steps', \
     file index 7 using 1:3 w lp title '2200 Steps', \
     file index 8 using 1:3 w lp title '14500 Steps'