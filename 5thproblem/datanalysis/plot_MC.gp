set terminal png size 750,462
set output "RDF_MC.png"


plot "out.5000" u 1:2 w lp title "5000 steps",\
 "out.10000" u 1:2 w lp title "10000 steps",\
	"out.100000" u 1:2 w lp title "100000 steps",\
	"out.1000000" u 1:2 w lp title "1000000 steps"

