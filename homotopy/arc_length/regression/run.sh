#!/bin/bash
debug="all"
#./homotopy_solver -p ../test_functions/f1.so -t none -d $debug -e jacobian -m 10 > f1_none
#./homotopy_solver -p ../test_functions/f1.so -t difference -d $debug -e jacobian > f1_difference
#./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -m 10 > f1_differential
./homotopy_solver -p ../test_functions/f1.so -t differential -d $debug -e jacobian -o f1.raw -m 100 -k exact -y 100 -w 0.05 -q tangent > f1_differential_exact
./homotopy_solver -p ../test_functions/f2.so -t differential -d $debug -e jacobian -o f2.raw -m 100 -k exact -y 100 -w 0.05 -q tangent > f2_differential_exact
./homotopy_solver -p ../test_functions/f3.so -t differential -d $debug -e jacobian -o f3.raw -m 100 -k exact -y 1000 -w 0.25 -q tangent > f3_differential_exact
./homotopy_solver -p ../test_functions/f3.so -t differential -d $debug -e jacobian -o f3.fail.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f3_differential_exact
./homotopy_solver -p ../test_functions/f3.so -t differential -d $debug -e jacobian -o f3.finegrid.raw -m 100 -k exact -y 1000 -w 0.01 -q tangent > f3_differential_exact
./homotopy_solver -p ../test_functions/f4.so -t differential -d $debug -e jacobian -o f4.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f4_differential_exact
./homotopy_solver -p ../test_functions/f5.so -t differential -d $debug -e jacobian -o f5.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f5_differential_exact
./homotopy_solver -p ../test_functions/f5.so -t differential -d $debug -e jacobian -o f5.sub.raw -m 100 -s sub_cross_product -k exact -y 1000 -w 0.20 -q tangent > f5_differential_exact
./homotopy_solver -p ../test_functions/f5.so -t differential -d $debug -e jacobian -o f5.diag.raw -m 100 -s diag_cross_product -k exact -y 1000 -w 0.20 -q tangent > f5_differential_exact
./homotopy_solver -p ../test_functions/f6.so -t differential -d $debug -e jacobian -o f6.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f6.so -t differential -d $debug -e jacobian -o f6.sub.raw -m 100 -s sub_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f6.so -t differential -d $debug -e jacobian -o f6.diag.raw -m 100 -s diag_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f7.so -t differential -d $debug -e jacobian -o f7.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f7.so -t differential -d $debug -e jacobian -o f7.sub.raw -m 100 -s sub_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f7.so -t differential -d $debug -e jacobian -o f7.diag.raw -m 100 -s diag_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f8.so -t differential -d $debug -e jacobian -o f8.raw -m 100 -k exact -y 1000 -w 0.05 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f8.so -t differential -d $debug -e jacobian -o f8.sub.raw -m 100 -s sub_cross_product -k exact -y 1000 -w 0.05 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f8.so -t differential -d $debug -e jacobian -o f8.diag.raw -m 100 -s diag_cross_product -k exact -y 1000 -w 0.05 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f9.so -t differential -d $debug -e jacobian -o f9.raw -m 100 -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f9.so -t differential -d $debug -e jacobian -o f9.sub.raw -m 100 -s sub_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
./homotopy_solver -p ../test_functions/f9.so -t differential -d $debug -e jacobian -o f9.diag.raw -m 100 -s diag_cross_product -k exact -y 1000 -w 0.20 -q tangent > f6_differential_exact
