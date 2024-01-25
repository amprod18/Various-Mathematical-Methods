set terminal png size 640, 480; set output "figE1.png"

set title "f(theta) = 3 - |theta - 4·cos(theeta)|"
set xlabel "theta (rad)"
set ylabel "f(theta)"
set grid


file = "res1a.dat"

plot file index 0 using 1:2 w lp notitle