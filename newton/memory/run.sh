#!/bin/bash

dynamic_analyis="valgrind --leak-check=full --show-leak-kinds=all -v"

all_case="f1 f2 f3 f4"
for case in ${all_case}; do
	${dynamic_analyis} ./newton_solver -i normal -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_jacobian -p ../test_functions/${case}.so >& ${case}_jacobian.valgrind
	${dynamic_analyis} ./newton_solver -i normal -e forward -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_forward -p ../test_functions/${case}.so >& ${case}_forward.valgrind
	${dynamic_analyis} ./newton_solver -i normal -e central -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_central -p ../test_functions/${case}.so >& ${case}_central.valgrind
	${dynamic_analyis} ./newton_solver -i chord -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_chord -p ../test_functions/${case}.so >& ${case}_chord.valgrind
	${dynamic_analyis} ./newton_solver -i broyden -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden -p ../test_functions/${case}.so >& ${case}_broyden.valgrind
	${dynamic_analyis} ./newton_solver -i broyden_inverted -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted -p ../test_functions/${case}.so >& ${case}_broyden_inverted.valgrind
	${dynamic_analyis} ./newton_solver -i broyden_inverted_bad -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted_bad -p ../test_functions/${case}.so >& ${case}_broyden_inverted_bad.valgrind
done

