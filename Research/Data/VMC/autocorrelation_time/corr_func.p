#set terminal epslatex
#set output 'optimizers.tex'
set terminal pdf
set output 'corr_func.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Correlation function"
set xlabel "k"
set ylabel "C(k)"
set xr [0:40]
set yr [0:1]
#set key off
plot    "./h20/corr_func.txt" using 1:2 title 'H_{20}' with linespoints, \
        "./h40/corr_func.txt" using 1:2 title 'H_{40}' with linespoints, \
        "./h60/corr_func.txt" using 1:2 title 'H_{60}' with linespoints, \
        "./h80/corr_func.txt" using 1:2 title 'H_{80}' with linespoints, \
        "./h100/corr_func.txt" using 1:2 title 'H_{100}' with linespoints, \
        "./h120/corr_func.txt" using 1:2 title 'H_{120}' with linespoints
