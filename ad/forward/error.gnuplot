# plot derivative error
#set mxtics 2
#set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set logscale x 10
set logscale y 10
set xrange [1:1e-20] 
set yrange [1e-20:1]  
set format x "10^{%T}"
set format y "10^{%T}"

set style line 1 lc "#000000" lw 5
set style line 2 lc "#FF0000" lw 3
set style line 3 lc "#00FF00" lw 3
set style line 4 lc "#0000FF" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced

epsilon=2.2204460493e-16

set output 'error1.png'
plot \
epsilon title 'epsilon' with lines linestyle 1,\
'data' using ($1):($2) title 'fd' with lines linestyle 2,\
'data' using ($1):($3) title 'cd' with lines linestyle 3,\
'data' using ($1):($4) title 'ad' with lines linestyle 4

set output 'error2.png'
plot \
epsilon title 'epsilon' with lines linestyle 1,\
'data' using ($1):($5) title 'fd2' with lines linestyle 3,\
'data' using ($1):($6) title 'ad2' with lines linestyle 4
