#gnuplot
set term pngcairo size 9cm,9.6cm font ",12"
set key autotitle columnhead
set output 'N2_py_sr_zoom.png'

unset key
set format y "%4.2f"
set format x "%3.1f"
set xlabel "R(a.u)"
#set ylabel "E(a.u)"
set ytics .02
set xtics .1
set border lw 3
plot [1.00:1.20][-0.01:+0.068]\
     'res_fci'    u 1:($2--109.278339)        w l dt 1 lc rgb 'gray'         lw 4 notitle,\
     'res_py'     u 1:($2--109.278339)        w l dt 1 lc rgb 'black'        lw 2 ,\
     'res_py'     u 1:($6+$2--109.278339)     w l dt 1 lc rgb 'light-red'    lw 2 ,\
     'res_py'     u 1:($7+$2--109.278339)     w l dt 1 lc rgb 'blue'         lw 2 ,\
     'res_py'     u 1:($7+$8+$2--109.278339)  w l dt 3 lc rgb 'dark-blue'    lw 2 ,\
     'res_molpro' u 1:($3--109.278339)        w l dt 1 lc rgb 'dark-green'   lw 2
