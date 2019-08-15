set terminal epslatex
set output 'corr.tex'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Autocorrelation time vs. Number of electrons"
set xlabel "Number of electrons"
set ylabel "Autocorrelation time"
set xr [10.0:110.0]
set yr [5:15]
set key off
plot    "data.txt" using 1:2 with linespoints
