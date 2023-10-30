set terminal png
set output "Energies.png"
file="SIM-L-032-TEMP-1500-MCTOT-10K.dat"
file1="SIM-L-032-TEMP-1800-MCTOT-10K.dat"
file2="SIM-L-032-TEMP-2500-MCTOT-10K.dat"
file3="SIM-L-032-TEMP-3500-MCTOT-10K.dat"
file4="SIM-L-032-TEMP-4500-MCTOT-10K.dat"

set ylabel "Energia"
set xlabel "N"
set grid
set key outside right top box

set style data lines

plot file using 1:2 title "T=1.5", \
file1 using 1:2 title "T=1.8", \
file2 using 1:2 title "T=2.5", \
file3 using 1:2 title "T=3.5", \
file4 using 1:2 title "T=4.5"
set output
