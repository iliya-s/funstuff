#gnuplot
set term pngcairo size 12.8cm,9.6cm font ",12"
set key autotitle columnhead
set output 'N2_py_mrpt.png'

set key left font ",10" samplen 2 Left reverse opaque width 0
set format y "%3.1f"
set format x "%3.1f"
set xlabel "R(a.u)"
set ylabel "E(a.u)"
set ytics .4
set xtics .5
set border lw 3
plot [0.5:3][-0.12:1.27]\
     'res_fci'    u 1:($2--109.278339)       w l dt 1 lc rgb 'gray'         lw 4 notitle,\
     'res_py'     u 1:($9--109.278339)       w l dt 1 lc rgb 'black'        lw 2 ,\
     'res_py'     u 1:($9+$11--109.278339)   w l dt 1 lc  6                 lw 2
