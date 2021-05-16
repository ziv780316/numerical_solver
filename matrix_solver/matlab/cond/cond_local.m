function [cond_num, A_norm, A_inv_norm_local, optima_arg, x0, df_dx] = cond_local( A, max_explored_set, max_explored_iter, debug )
% cond_local
% use gradient search for condition number norm 1 evaluation

n = size(A, 1);

if nargin < 2
	max_explored_set = 1;
end
if nargin < 3
	max_explored_iter = n;
end
if nargin < 4
	debug = 0;
end

format long e;

A_norm = norm(A, 1);
A_inv = inv(A);
A_inv_norm_exact = norm(A_inv, 1);

argmax = 0;
f = 0;
optima_arg = 0;
optima = 0;
explored_set = zeros(n, 1);
x0 = zeros(n, 1);
total_iter = 0;
total_explored_set = 0;

for k = 1:1:max_explored_set
	n_unexplored = n - sum(explored_set);
	x0_old = x0;
	x0 = (1/(n_unexplored)) * (ones(n, 1) - explored_set);
	if x0 == x0_old
		fprintf( 'gradient search converged, x0 is the same\n' );
		break
	end
	x = x0;
	total_explored_set = total_explored_set + 1;
	for i = 1:1:max_explored_iter
		total_iter = total_iter + 1;
		fold = f;
		a_inv = A \ x;
		f = norm(a_inv, 1);
		if f > optima
			optima_arg = argmax;
			optima = f;
		end

		df_dg = sign(a_inv);
		df_dx = A' \ df_dg;
		argmax_old = argmax;
		%[df_dx_max, argmax] = max(abs(df_dx));
		[df_dx_max, argmax] = max(df_dx);
		fprintf( 'f%d%d=%.10e df_dx_max=%.10e argmax=%d f_optima=%.10e argmax_optima=%d\n', k, i, f, df_dx_max, argmax, optima, optima_arg );
		x = zeros(n, 1);
		x(argmax) = 1;
		explored_set(argmax) = 1;

		if i > 1
			if argmax == argmax_old
				fprintf( 'gradient search converged, argmax is the same\n' );
				break
			end
			%if f < f_old
			%	fprintf( 'f=%.10e < f_old=%.10e, stop\n', f, f_old );
			%	f = f_old;
			%	argmax = argmax_old;
			%	break
			%end
		end
	end
end

fprintf( 'spend %d iter to explored %d set, local max f=%.10e argmax=%d\n', total_iter, total_explored_set, optima, optima_arg );
cond_num = A_norm * optima;
A_inv_norm_local = optima;
fprintf( 'A_inv_norm_exact=%.10e A_inv_norm_local=%.10e\n', A_inv_norm_exact, optima );
fprintf( 'cond_exact=%.10e cond_local=%.10e\n', A_norm * A_inv_norm_exact, A_norm * optima );

end
