set terminal epslatex lw 5
set output 'eloc.tex'
#set terminal pdf lw 2
#set output 'eloc.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "$KS_z$GHF"
set xlabel "$theta$"
set ylabel "Energy, Ha"
set xr [0:3.1415]
set yr [-10.0:10.0]
set key top left
plot    "ksghf.txt" using 1:($3) title 'E_L' with linespoints pointtype -1, \
        "ksghf.txt" using 1:($4) title 'T_L' with linespoints pointtype -1, \
        "ksghf.txt" using 1:($5) title 'V_L' with linespoints pointtype -1,
