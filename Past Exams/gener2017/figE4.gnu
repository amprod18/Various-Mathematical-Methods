set terminal png size 640, 480; set output "figE4.png"

set title "Integral2 Convergence"
set xlabel "N Samples"
set ylabel "|I2_{computed} - I2_{real}|"
set grid
set logscale x

file = "res1a.dat"

plot file index 4 using 1:4 w lp notitle