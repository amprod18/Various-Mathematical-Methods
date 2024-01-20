set terminal png size 640, 480; set output 'Root_Methods_fig4.png'
set autoscale
set ylabel "Abnormal Excentricity"
set xlabel "t (years)"
set grid
set title "Abnormal Excentricity of comet Hale-Bopp"

file = 'Root_Methods_res.dat'

plot file using 1:2 index 2 title 'Newton-Raphson',\
     file using 1:5 index 2 title 'Secant'