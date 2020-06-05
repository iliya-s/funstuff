set terminal epslatex lw 5
set output 'cusp.tex'
#set terminal pdf lw 2
#set output 'cusp.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Electron Cusp"
set xlabel "$theta$"
set ylabel "$Psi$"
set xr [0:3.1415]
set yr [:]
set key right bottom
plot    "ksghf.txt" using 1:($2/0.192224) title 'KSGHF' with linespoints pointtype -1, \
        "jksghf.txt" using 1:($2/0.298160) title 'JKSGHF' with linespoints pointtype -1,
