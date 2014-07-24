reset
set term pdf font "Times,2"    #set terminal and output file
set output "distribution_taille_seq.pdf"
set xlabel "x value"    #set x and y label
set ylabel "Frequency"
set xrange [0:14000]    #set x and y range
set yrange [0:1400000]
set style fill solid    #set plot style	
set boxwidth 0.5
unset key    #legend not plotted
plot "distribution_taille_seq.dat" u 1:2 w linespoints #using 2:xtic(1) with boxes
#plot bar chart and the value labels on the bars
