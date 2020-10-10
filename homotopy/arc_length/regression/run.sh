#!/bin/bash
debug="all"
#./homotopy_solver -p ../test_functions/f1.so -t none -d $debug -e jacobian -m 10 > f1_none
#./homotopy_solver -p ../test_functions/f1.so -t difference -d $debug -e jacobian > f1_difference
#./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -m 10 > f1_differential
./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -o f1.raw -m 100 -k exact -y 100 -w 0.2 -q tangent > f1_differential_exact
