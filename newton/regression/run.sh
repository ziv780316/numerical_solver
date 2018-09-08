#!/bin/bash

# basic test
all_case="f1 f2 f3 f4"
for case in ${all_case}; do
	./newton_solver -i normal -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_jacobian -p ../test_functions/${case}.so
	./newton_solver -i normal -e forward -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_forward -p ../test_functions/${case}.so
	./newton_solver -i normal -e central -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_central -p ../test_functions/${case}.so
	./newton_solver -i chord -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_chord -p ../test_functions/${case}.so
	./newton_solver -i broyden -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden -p ../test_functions/${case}.so
	./newton_solver -i broyden_inverted -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted -p ../test_functions/${case}.so
	./newton_solver -i broyden_inverted_bad -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted_bad -p ../test_functions/${case}.so
done

# modified newton in ill-conditioned problem
./newton_solver -i normal -e jacobian -j 1e-12 -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f5_jacobian -p ../test_functions/f5.so
./newton_solver -i normal -e jacobian -f damped -j 1e-12 -b 1.0 -x f2_x0_1 -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f2_jacobian_damped -p ../test_functions/f2.so
./newton_solver -i normal -e jacobian -j 0     -n 10 -b 1.0 -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f6_jacobian -p ../test_functions/f6.so
./newton_solver -i normal -e jacobian -j 1e-12 -n 10 -b 1.0 -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f6_jacobian_jmin -p ../test_functions/f6.so
./newton_solver -i normal -e jacobian -j 1e-3  -n 10 -b 1.0 -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f6_jacobian_jmin_3 -p ../test_functions/f6.so

# compare with golden
for file in `ls golden`; do
	diff -q ${file/golden\//} $file
done
