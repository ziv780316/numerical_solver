
function x = qr_linsolve( A, b, debug )

	size_A = size(A);
	nrow = size_A(1);
	ncol = size_A(2);

	[Q,R,P,rank]=qr_householder(A,0,0);

	R11 = R(1:rank,1:rank);
	y = Q'*b;
	y1 = y(1:rank);
	x1 = R11\y1;
	x = zeros(ncol,1);
	x(1:rank) = x1;
	x = P*x;

	if debug
		fprintf('Q:\n');
		disp(Q);
		fprintf('R11:\n');
		disp(R11);
		fprintf('P:\n');
		disp(P);
		fprintf('y1:\n');
		disp(y1);
		fprintf('x1:\n');
		disp(x1);
		fprintf('x:\n');
		disp(x);
		fprintf('A*x:\n');
		disp(A*x);
	end

	%C = R' * R;
	%y = C\b;
	%x = Q * R * y;
end
