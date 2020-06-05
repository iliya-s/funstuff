set terminal epslatex lw 5
set output 'jeloc.tex'
#set terminal pdf lw 2
#set output 'jeloc.pdf'
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "J-$KS_z$GHF"
set xlabel "$theta$"
set ylabel "Energy, Ha"
set xr [0:3.1415]
set yr [-10.0:10.0]
set key top left
plot    "jksghf.txt" using 1:($3) title 'E_L' with linespoints pointtype -1, \
        "jksghf.txt" using 1:($4) title 'T_L' with linespoints pointtype -1, \
        "jksghf.txt" using 1:($5) title 'V_L' with linespoints pointtype -1,
