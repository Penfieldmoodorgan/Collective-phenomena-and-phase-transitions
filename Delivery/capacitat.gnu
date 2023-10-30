set terminal png
set output "Capacitat.png"
file="MC-L-064-DADES-CAPACITAT.dat"


set ylabel "Cv"
set xlabel "T"
set grid
set key outside right top box

set style data lines

plot file using 1:2 title "Cv"
set output
