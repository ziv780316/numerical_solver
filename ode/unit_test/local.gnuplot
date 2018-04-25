
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:25] 
set yrange [0:1.5] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3
set style line 7 lc "#FF3399" lw 3
set style line 8 lc "#664443" lw 3
set style line 9 lc "#004443" lw 3
set style line 10 lc "#6644FF" lw 3

set term png size 960,840 font 20 enhanced
set output 'local.png'

step(t,a) = (t >= a) ? 1 : 1/0
f(t)=exp(t)/(99+exp(t))
f0(t)=((1.0102040816e-02)*exp(9.8000000000e-01*(t - 0.0000000000e+00)) - 1.0102040816e-02 + 1.0000000000e-02)*step(t,0.0000000000e+00)
#f2500(t)=((1.2480036646e-01)*exp(7.8107890977e-01*(t - 2.5000000000e+00)) - 1.2480036646e-01 + 1.0946054511e-01)*step(t,2.5000000000e+00)
#f5000(t)=((-1.2068346736e+00)*exp(-1.9895380159e-01*(t - 5.0000000000e+00)) - -1.2068346736e+00 + 5.9947690080e-01)*step(t,5.0000000000e+00)
#f7500(t)=((-5.4964343436e-02)*exp(-8.9609532686e-01*(t - 7.5000000000e+00)) - -5.4964343436e-02 + 9.4804766343e-01)*step(t,7.5000000000e+00)
#f10000(t)=((-4.4928869186e-03)*exp(-9.9105459741e-01*(t - 1.0000000000e+01)) - -4.4928869186e-03 + 9.9552729871e-01)*step(t,1.0000000000e+01)
f2(t)=((4.1089495730e-02)*exp(9.2119202000e-01*(t - 2.0000000000e+00)) - 4.1089495730e-02 + 3.9403990000e-02)*step(t,2.0000000000e+00)
f4(t)=((1.7993261999e-01)*exp(7.0291554219e-01*(t - 4.0000000000e+00)) - 1.7993261999e-01 + 1.4854222891e-01)*step(t,4.0000000000e+00)
f6(t)=((4.8706842995e+00)*exp(5.1192975051e-02*(t - 6.0000000000e+00)) - 4.8706842995e+00 + 4.7440351247e-01)*step(t,6.0000000000e+00)
f8(t)=((-8.3187986901e-02)*exp(-8.4737003219e-01*(t - 8.0000000000e+00)) - -8.3187986901e-02 + 9.2368501609e-01)*step(t,8.0000000000e+00)
f10(t)=((-3.3919855959e-05)*exp(-9.9993216259e-01*(t - 1.0000000000e+01)) - -3.3919855959e-05 + 9.9996608129e-01)*step(t,1.0000000000e+01)
f12(t)=((-0.0000000000e+00)*exp(-1.0000000000e+00*(t - 1.2000000000e+01)) - -0.0000000000e+00 + 1.0000000000e+00)*step(t,1.2000000000e+01)

plot f(x) title 'exact' with lines linestyle 1,\
'output_ivp4_ab1_t1p0' using ($1):($2)  title 'fe' with lines linestyle 10,\
f0(x) title 'f0' with lines linestyle 2,\
f2(x) title 'f2' with lines linestyle 3,\
f4(x) title 'f4' with lines linestyle 4,\
f6(x) title 'f6' with lines linestyle 6,\
f8(x) title 'f8' with lines linestyle 7,\
f10(x) title 'f10' with lines linestyle 8,\
f12(x) title 'f12' with lines linestyle 9
