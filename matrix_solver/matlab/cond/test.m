A = readmatrix( 'ill_A' );
max_explored_set = 5;
max_explored_iter = size(A, 1);
debug = 0;
[cond_num, A_norm, A_inv_norm_local, optima_arg, x0, df_dx] = cond_local( A, max_explored_set, max_explored_iter, debug );
