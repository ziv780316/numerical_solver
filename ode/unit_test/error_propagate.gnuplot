
set mxtics 2
set mytics 2
set grid
set grid mxtics
set grid mytics
set grid xtics ytics
set xrange [0:10] 
#set yrange [-1.5:1.5] 

set style line 1 lc "#000000" lw 3
set style line 2 lc "#0000FF" lw 3
set style line 3 lc "#FF0000" lw 3
set style line 4 lc "#00FF00" lw 3
set style line 5 lc "#FF00FF" lw 3
set style line 6 lc "#FF8833" lw 3

set term png size 960,840 font 20 enhanced
set output 'error_propagate.png'

l = 1
h = 0.1
d = 0.05
step(x,a) = x>a ? 1 : 0
f(x) = exp(-l*x)
f1(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h))
f2(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h))
f3(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h))
f4(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h))
f5(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h))
f6(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h))
f7(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h))
f8(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h))
f9(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h))
f10(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h))
f11(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h))
f12(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h))
f13(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h))
f14(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h))
f15(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h))
f16(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h)) - (d * exp(-l*(x-16*h)) * step(x,16*h))
f17(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h)) - (d * exp(-l*(x-16*h)) * step(x,16*h)) - (d * exp(-l*(x-17*h)) * step(x,17*h))
f18(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h)) - (d * exp(-l*(x-16*h)) * step(x,16*h)) - (d * exp(-l*(x-17*h)) * step(x,17*h)) - (d * exp(-l*(x-18*h)) * step(x,18*h))
f19(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h)) - (d * exp(-l*(x-16*h)) * step(x,16*h)) - (d * exp(-l*(x-17*h)) * step(x,17*h)) - (d * exp(-l*(x-18*h)) * step(x,18*h)) - (d * exp(-l*(x-19*h)) * step(x,19*h))
f20(x) = exp(-l*x) - (d * exp(-l*(x-1*h)) * step(x,1*h)) - (d * exp(-l*(x-2*h)) * step(x,2*h)) - (d * exp(-l*(x-3*h)) * step(x,3*h)) - (d * exp(-l*(x-4*h)) * step(x,4*h)) - (d * exp(-l*(x-5*h)) * step(x,5*h)) - (d * exp(-l*(x-6*h)) * step(x,6*h)) - (d * exp(-l*(x-7*h)) * step(x,7*h)) - (d * exp(-l*(x-8*h)) * step(x,8*h)) - (d * exp(-l*(x-9*h)) * step(x,9*h)) - (d * exp(-l*(x-10*h)) * step(x,10*h)) - (d * exp(-l*(x-11*h)) * step(x,11*h)) - (d * exp(-l*(x-12*h)) * step(x,12*h)) - (d * exp(-l*(x-13*h)) * step(x,13*h)) - (d * exp(-l*(x-14*h)) * step(x,14*h)) - (d * exp(-l*(x-15*h)) * step(x,15*h)) - (d * exp(-l*(x-16*h)) * step(x,16*h)) - (d * exp(-l*(x-17*h)) * step(x,17*h)) - (d * exp(-l*(x-18*h)) * step(x,18*h)) - (d * exp(-l*(x-19*h)) * step(x,19*h)) - (d * exp(-l*(x-20*h)) * step(x,20*h))


plot f(x) title 'f' with lines linestyle 2,\
f1(x) title 'f1' with lines linestyle 1,\
f2(x) title 'f2' with lines linestyle 1,\
f3(x) title 'f3' with lines linestyle 1,\
f4(x) title 'f4' with lines linestyle 1,\
f5(x) title 'f5' with lines linestyle 1,\
f6(x) title 'f6' with lines linestyle 1,\
f7(x) title 'f7' with lines linestyle 1,\
f8(x) title 'f8' with lines linestyle 1,\
f9(x) title 'f9' with lines linestyle 1,\
f10(x) title 'f10' with lines linestyle 1,\
f11(x) title 'f11' with lines linestyle 1,\
f12(x) title 'f12' with lines linestyle 1,\
f13(x) title 'f13' with lines linestyle 1,\
f14(x) title 'f14' with lines linestyle 1,\
f15(x) title 'f15' with lines linestyle 1,\
f16(x) title 'f16' with lines linestyle 1,\
f17(x) title 'f17' with lines linestyle 1,\
f18(x) title 'f18' with lines linestyle 1,\
f19(x) title 'f19' with lines linestyle 1,\
f20(x) title 'f20' with lines linestyle 1
