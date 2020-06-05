set terminal epslatex linewidth 2
set output 'C2.tex'
#set terminal pdf linewidth 2
#set output 'C2.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Carbon Dimer Optimization"
set xlabel "Iteration"
set ylabel "Energy, Ha"
set xr [0:40]
set yr [:]
set key right
plot    "data" using 0:1 title 'aLM-1' with linespoints, \
        "data" using 0:2 title 'aLM-2' with linespoints, \
        "data" using 0:3 title 'aLM-3' with linespoints, \
        "data" using 0:4 title 'LM' with linespoints, \
        "data" using 0:5 title 'AMSGrad' with linespoints, \
        "data" using 0:6 title 'SR' with linespoints,
