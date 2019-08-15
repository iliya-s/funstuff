set terminal epslatex
set output 'scaling.tex'
#set terminal pdf
#set output 'scaling.pdf'
set autoscale                        # scale axes automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set log xy
#set title ""
set xlabel "Log of number of electrons"
set ylabel "Log of CPU time"
#set xr [18:180]
#set yr [1:100000]
#set key off
set key right center spacing 2
f(x) = a * x + b
fit f(x) 'hchain.dat' u (log($1)):(log($2)) via a, b
f(x) = a1 * x + b1
fit f(x) 'hchain.dat' u (log($1)):(log($3)) via a1,b1
set label 1 sprintf("slope = %3.4f", a) at 46,1 rotate by 9
set label 2 sprintf("slope = %3.4f", a1) at 46,740 rotate by 23
plot    "hchain.dat" using 1:3 title '$H_{x}$ old' with linespoints, \
        "hchain.dat" using 1:2 title '$H_{x}$ new' with linespoints
