set terminal epslatex
set output "N2curve.tex"
#set terminal pdf
#set output "N2curve.pdf"
set autoscale
set xtic auto
set ytic auto
set title "N_{2} Dissociation curve"
set xlabel "Bond length, Angstrom"
set ylabel "Energy, E_h"
set key bottom right
plot [0.7:2.7][-109.5:-108.5]\
     'res_fci'    using 1:($2)        title "Exact Diagonalization" with lines,\
     'res_py'     using 1:($2)        title "Hartree Fock" with lines,\
     'res_py'     using 1:($7+$2)     title "Coupled Cluster Singles and Doubles" with lines,\
     'res_molpro' using 1:($3)        title "Second Order Pertubation Theory" with lines
