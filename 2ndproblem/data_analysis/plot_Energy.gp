set terminal png size 750,462
set output "Energy_vs_step.png"

set yrange [-4.54:-4.30]


set xlabel "time"
set	ylabel "Total Energy" 

plot "ETP.out" u ($1*0.0005):6 w lp title "seed 1928427217",\
 "out.2nd" u ($1*0.0005):6 w lp title "seed 1928427307",\
 "out.3" u ($1*0.0005):6 w lp title "seed 1928427547"
