set terminal png
set output "Magnetitzacio.png"
file="MC-L-064-DADES-MAG.dat"


set ylabel "M/N"
set xlabel "T"
set grid
set key outside right top box

set style data lines

plot file using 1:2:4 with yerrorbars title "<|m|>", \
file using 1:3 title "sqrt(<m^2>)"
set output
