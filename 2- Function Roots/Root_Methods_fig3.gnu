set terminal png size 640, 480; set output 'Root_Methods_fig3.png'
set autoscale
set ylabel "Orbital Components"
set xlabel "E (rad)"
set grid
set title "Distances of comet Hale-Bopp"
set key

file = 'Root_Methods_res.dat'

plot file using 2:3 index 2 title 'x NR',\
     file using 2:4 index 2 title 'y NR',\
     file using 2:6 index 2 title 'x S',\
     file using 2:7 index 2 title 'y S'