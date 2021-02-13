
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

set yrange [-1.5:1] 
set output 'f3.png'
plot 'f3.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f3.fail.png'
plot 'f3.fail.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f3.finegrid.png'
plot 'f3.finegrid.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f4.png'
plot 'f4.raw.trace' using ($1):($2) title 'x' with linespoints linestyle 2 pointtype 7 pointsize 2

set output 'f5.png'
plot 'f5.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f5.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f5.sub.png'
plot 'f5.sub.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f5.sub.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f5.diag.png'
plot 'f5.diag.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f5.diag.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f6.png'
plot 'f6.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f6.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f6.sub.png'
plot 'f6.sub.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f6.sub.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f6.diag.png'
plot 'f6.diag.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f6.diag.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f7.png'
plot 'f7.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f7.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f7.sub.png'
plot 'f7.sub.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f7.sub.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f7.diag.png'
plot 'f7.diag.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f7.diag.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f8.png'
plot 'f8.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f8.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f8.sub.png'
plot 'f8.sub.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f8.sub.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f8.diag.png'
plot 'f8.diag.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f8.diag.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2

set output 'f9.png'
plot 'f9.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f9.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2 ,\
'f9.raw.trace' using ($1):($4) title 'x3' with linespoints linestyle 4 pointtype 7 pointsize 2

set output 'f9.sub.png'
plot 'f9.sub.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f9.sub.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2 ,\
'f9.sub.raw.trace' using ($1):($4) title 'x3' with linespoints linestyle 4 pointtype 7 pointsize 2

set output 'f9.diag.png'
plot 'f9.diag.raw.trace' using ($1):($2) title 'x1' with linespoints linestyle 2 pointtype 7 pointsize 2, \
'f9.diag.raw.trace' using ($1):($3) title 'x2' with linespoints linestyle 3 pointtype 7 pointsize 2, \
'f9.diag.raw.trace' using ($1):($4) title 'x3' with linespoints linestyle 4 pointtype 7 pointsize 2
