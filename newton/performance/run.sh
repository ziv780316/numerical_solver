#!/bin/bash

bench_llc=0
perf_sw_evnet="cpu-clock,context-switches,cpu-migrations"
if [ $bench_llc -eq 1 ]; then
	perf_hw_event="instructions,cache-misses,cache-references,LLC-load-misses,LLC-loads"
else
	perf_hw_event="instructions,cache-misses,cache-references,L1-dcache-load-misses,L1-dcache-loads"
fi
bench="perf stat --event=${perf_hw_event},${perf_sw_evnet}"

all_case="f1 f2 f3 f4"
for case in ${all_case}; do
	${bench} ./newton_solver -i normal -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_jacobian -p ../test_functions/${case}.so >& ${case}_jacobian.perf
	${bench} ./newton_solver -i normal -e forward -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_forward -p ../test_functions/${case}.so >& ${case}_forward.perf
	${bench} ./newton_solver -i normal -e central -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_central -p ../test_functions/${case}.so >& ${case}_central.perf
	${bench} ./newton_solver -i chord -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_chord -p ../test_functions/${case}.so >& ${case}_chord.perf
	${bench} ./newton_solver -i broyden -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden -p ../test_functions/${case}.so >& ${case}_broyden.perf
	${bench} ./newton_solver -i broyden_inverted -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted -p ../test_functions/${case}.so >& ${case}_broyden_inverted.perf
	${bench} ./newton_solver -i broyden_inverted_bad -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o ${case}_broyden_inverted_bad -p ../test_functions/${case}.so >& ${case}_broyden_inverted_bad.perf
done

