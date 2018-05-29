
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
set style line 7 lc "#FF3399" lw 3
set style line 8 lc "#664443" lw 3

set term png size 960,840 font 20 enhanced
set output 'rk_ivp4.png'

f(x)=exp(x)/(99+exp(x))

plot f(x) title 'exact' with lines linestyle 1,\
'output_ivp4_rk4_t0p1' title 'RK4\_0p1' with lines linestyle 2

