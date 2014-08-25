reset
set term pdf font "Times,5"    #set terminal and output file
set output "true_positive.pdf"
set xlabel "read length"    #set x and y label
set ylabel "true positive"
set xrange [0:300]    #set x and y range
set yrange [0:50]
set style fill solid    #set plot style	
set boxwidth 0.5
#unset key    #legend not plotted
set arrow 1 from 0,42 to 300,42 
plot "result.dat" u 1:2 smooth bezier title 'complete_genome'
plot  42   with lines lt 1 title "y = 42"
