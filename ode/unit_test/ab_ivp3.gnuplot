
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:25] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced
set output 'ab_ivp3.png'

f(x)=1/(x+1)

plot f(x) title '1/(x+1)' with lines linestyle 1,\
'output_ivp3_ab1_t0p19' using ($1):($2) title 'AB1_0p19' with lines linestyle 2,\
'output_ivp3_ab1_t0p21' using ($1):($2) title 'AB1_0p21' with lines linestyle 3


