set terminal png
set output "Susceptibilitat.png"
file="MC-L-064-DADES-SUSCEPTIBILITA.dat"


set ylabel "X"
set xlabel "T"
set grid
set key outside right top box

set style data lines

plot file using 1:2 title "<X*>"
set output
