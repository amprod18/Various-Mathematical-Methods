set terminal png size 640, 480; set output 'Root_Methods_fig1.png'
set autoscale
set ylabel "Distance(E) / Distance'(E)"
set xlabel "E (rad)"
set grid
set title "Distance and velocity of comet Hale-Bopp"
set key left top

file = 'Root_Methods_res.dat'

plot file using 1:2 index 0 title "Distance(E)", \
     file using 1:3 index 0 title "Distance'(E)"