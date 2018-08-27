#!/bin/bash

all_case="f1 f2 f3"
for case in ${all_case}; do
	./newton_solver -i normal -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_jacobian -p ../test_functions/${case}.so
	./newton_solver -i normal -e forward -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_forward -p ../test_functions/${case}.so
	./newton_solver -i normal -e central -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_central -p ../test_functions/${case}.so
	./newton_solver -i chord -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_chord -p ../test_functions/${case}.so
	./newton_solver -i broyden -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden -p ../test_functions/${case}.so
	./newton_solver -i broyden_inverted -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted -p ../test_functions/${case}.so
	./newton_solver -i broyden_inverted_bad -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted_bad -p ../test_functions/${case}.so

	for file in `ls ${case}_*`; do
		diff -q $file golden/$file
	done
done

