
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:1] 
#set yrange [-5:5] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced

set output 'f1.png'
plot 'f1.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f2.png'
plot 'f2.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f3.png'
plot 'f3.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f3.fail.png'
plot 'f3.fail.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f3.finegrid.png'
plot 'f3.finegrid.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2
