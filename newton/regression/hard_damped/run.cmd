###QA_TAG newton
###QA_GOLDEN_FILE run.log
###QA_GOLDEN_FILE run.log.nr

# paramters
test_so=${DEVELOP_ROOT_DIR}/numerical_solver/newton/test_functions/f3.so 
nr_method=normal
jacobian_matrix=jacobian
output=run.log 
nr_maxiter=100
rtol=1e-3
atol=1e-6
ftol=1e-9
damp=damped
max_dx=3

# run
newton_solver \
-i $nr_method \
-e $jacobian_matrix \
-m $nr_maxiter \
-r $rtol \
-a $atol \
-u $ftol \
-o $output \
-p $test_so \
-f $damp \
-b $max_dx \
-d 

gnuplot plot.gnuplot
