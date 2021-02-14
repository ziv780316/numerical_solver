###QA_TAG homotopy arc_length
###QA_GOLDEN_FILE run.log.cross
###QA_GOLDEN_FILE run.log.cross.tangent
###QA_GOLDEN_FILE run.log.cross.trace
###QA_GOLDEN_FILE run.log.diag
###QA_GOLDEN_FILE run.log.diag.tangent
###QA_GOLDEN_FILE run.log.diag.trace

# paramters
test_so=${DEVELOP_ROOT_DIR}/numerical_solver/homotopy/arc_length/test_functions/f5.so 
extrapolate_type=differential
debug_type=all
jacobian_matrix=jacobian
output=run.log.cross 
nr_maxiter=100
dfdp_type=exact
maxstep=1000
length=0.05
constrain=tangent
backtrace_handle=cross_product

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
-q $constrain \
-s $backtrace_handle

output=run.log.diag 
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
-q $constrain \

gnuplot plot.gnuplot
