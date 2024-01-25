set terminal png size 640, 480; set output "figE2.png"

set title "Simpsons Convergence"
set xlabel "Step Size (rad)"
set ylabel "I_{computed}"
set grid
set logscale x

file = "res1a.dat"

plot file index 2 using 1:2 w lp notitle