#!/bin/bash
debug="all"
#./homotopy_solver -p ../test_functions/f1.so -t none -d $debug -e jacobian -m 10 > f1_none
#./homotopy_solver -p ../test_functions/f1.so -t difference -d $debug -e jacobian > f1_difference
#./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -m 10 > f1_differential
./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -m 100 -k exact -y 20 -w 0.1 -q tangent > f1_differential_exact
