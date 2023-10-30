set terminal png
set output "Energia.png"
file="MC-L-064-DADES-ENERGIA.dat"


set ylabel "E/N"
set xlabel "T"
set grid
set key outside right top box

set style data lines

plot file using 1:2:3 with yerrorbars title "<e>"
set output
