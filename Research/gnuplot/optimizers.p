#set terminal epslatex
#set output 'optimizers.tex'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Various optimizers vs. Iteration"
set xlabel "Iteration"
set ylabel "Energy"
set xr [0:40]
set yr [-10.6:-10.1]
#set key off
plot    "data.txt" using 0:1 title 'AMSGrad' with linespoints, \
        "data.txt" using 0:2 title 'SGD' with linespoints, \
        "data.txt" using 0:3 title 'SR' with linespoints,
