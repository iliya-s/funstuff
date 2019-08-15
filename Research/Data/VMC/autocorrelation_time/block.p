set terminal epslatex
set output 'optimizers.tex'
#set terminal pdf
#set output 'block.pdf'
set autoscale                        # scale axes automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set log x
set title "Blocking"
set xlabel "Block length"
set ylabel "r_x"
set xr [1:20000]
set yr [1:50]
#set key off
set key left
plot    "./h20/block.txt" using 1:2 title 'H_{20}' with linespoints, \
        "./h40/block.txt" using 1:2 title 'H_{40}' with linespoints, \
        "./h60/block.txt" using 1:2 title 'H_{60}' with linespoints, \
        "./h80/block.txt" using 1:2 title 'H_{80}' with linespoints, \
        "./h100/block.txt" using 1:2 title 'H_{100}' with linespoints, \
        "./h120/block.txt" using 1:2 title 'H_{120}' with linespoints
