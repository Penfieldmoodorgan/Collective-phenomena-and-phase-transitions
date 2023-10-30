set terminal png
set output "P1-configuration.png"
file="configuration.dat"
L=64
symbsize=1.2
set size square
set xrange [0.5:L+0.5]
set yrange [0.5:L+0.5]
plot file using 1:2 with points pt 5 ps symbsize t " "
set output


