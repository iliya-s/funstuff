#gnuplot
set term pngcairo size 12.8cm,9.6cm font ",12"
set key autotitle columnhead
set output 'BH_py_mrpt.png'

set key left font ",10" samplen 2 Left reverse opaque width 0
set format y "%3.1f"
set format x "%3.1f"
set xlabel "R(a.u)"
set ylabel "E(a.u)"
set ytics .1
set xtics .5
set border lw 3
plot [0.5:3][-0.05:0.40]\
     'res_fci'    u 1:($2--25.235431)   w l dt 1 lc rgb 'gray'         lw 4 notitle,\
     'res_py'     u 1:($7--25.235431)   w l dt 1 lc rgb 'black'        lw 2 ,\
     'res_py'     u 1:($9--25.235431)   w l dt 1 lc  6                 lw 2
