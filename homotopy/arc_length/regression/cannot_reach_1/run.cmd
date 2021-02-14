###QA_TAG homotopy arc_length
###QA_GOLDEN_FILE run.log
###QA_GOLDEN_FILE run.log.tangent
###QA_GOLDEN_FILE run.log.trace

# paramters
test_so=${DEVELOP_ROOT_DIR}/numerical_solver/homotopy/arc_length/test_functions/f2.so 
extrapolate_type=differential
debug_type=all
jacobian_matrix=jacobian
output=run.log 
nr_maxiter=100
dfdp_type=exact
maxstep=100
length=0.05
constrain=tangent

# run
homotopy_solver \
-p $test_so \
-t $extrapolate_type \
-d $debug_type \
-e $jacobian_matrix \
-o $output \
-t $extrapolate_type \
-d $debug_type \
-e jacobian \
-o $output \
-m $nr_maxiter \
-k $dfdp_type \
-y $maxstep \
-w $length \
-q $constrain

gnuplot plot.gnuplot
