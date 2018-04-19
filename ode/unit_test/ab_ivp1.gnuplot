
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:100] 
#set yrange [-1.5:1.5] 
set yrange [-5:5] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced
set output 'ab_ivp1.png'

f(x)=exp(-x)

plot f(x) title 'e^-t' with lines linestyle 1,\
'output_ivp1_ab1_t1p9' using ($1):($2) title 'AB1_1p9' with lines linestyle 2,\
'output_ivp1_ab1_t2p1' using ($1):($2) title 'AB1_2p1' with lines linestyle 3,\
'output_ivp1_ab3_t0p5' using ($1):($2) title 'AB3_0p5' with lines linestyle 4,\
'output_ivp1_ab3_t0p6' using ($1):($2) title 'AB3_0p6' with lines linestyle 5


