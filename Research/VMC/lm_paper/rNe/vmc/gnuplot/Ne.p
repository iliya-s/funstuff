set terminal epslatex linewidth 2
set output 'Ne.tex'
#set terminal pdf linewidth 2
#set output 'Ne.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Neon Atom Optimization"
set xlabel "Iteration"
set ylabel "Energy, Ha"
set xr [0:10]
set yr [:-128.8]
set key right
plot    "data" using 0:1 title 'aLM-1' with linespoints, \
        "data" using 0:2 title 'aLM-2' with linespoints, \
        "data" using 0:3 title 'aLM-3' with linespoints, \
        "data" using 0:4 title 'LM' with linespoints, \
