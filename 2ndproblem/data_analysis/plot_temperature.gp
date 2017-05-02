set terminal png size 750,462
set output "temp_vs_time.png"

#set yrange [4.10:0.60]


set xlabel "time"
set	ylabel "Temperature"

plot "ETP.out" u ($1*0.0005):2 w lp title "seed 1928427217",\
 "out.2nd" u ($1*0.0005):2 w lp title "seed 1928427307",\
 "out.3" u ($1*0.0005):2 w lp title "seed 1928427547"
