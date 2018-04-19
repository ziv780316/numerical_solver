
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:100] 
#set yrange [-1.5:1.5] 
set yrange [-1:1] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced
set output 'am.png'

f(x)=exp(-x)

plot f(x) title 'e^-t' with lines linestyle 1,\
'output_am1_t2p5' using ($1):($2) title 'AM1_2p5' with lines linestyle 2,\
'output_am2_t0p9' using ($1):($2) title 'AM2_0p9' with lines linestyle 3,\
'output_am2_t1p0' using ($1):($2) title 'AM2_1p0' with lines linestyle 4,\
'output_am3_t5p9' using ($1):($2) title 'AM3_5p9' with lines linestyle 5,\
'output_am3_t6p1' using ($1):($2) title 'AM3_6p1' with lines linestyle 6


