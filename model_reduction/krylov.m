rng(1);
m = 3;
A = rand(m,m) + eye(m,m);
b = rand(m,1);

D = diag(diag(A));
P = eye(m,m);
%P = D;
%P = full(ilu(sparse(A)));
A = P\A;
b = P\b;
x = A\b;
x_appr = 0;

for n=1:1:m
	[Q1, Q2, H] = arnoldi_mgs( A, b, n );
	y1 = (A*Q1)\b;
	y2 = H\(norm(b)*eye(n+1,1));
	y3 = (H'*H)\(H'*norm(b)*eye(n+1,1));
	x1 = Q1*y1;
	x2 = Q1*y2;
	x3 = Q1*y3;
	r1 = norm(b - A*x1);
	r2 = norm(b - A*x2);
	r3 = norm(b - A*x3);
	fprintf( 'r1=%e\n', r1 );
	fprintf( 'x1=\n' );
	disp(x1);
	fprintf( 'r2=%e\n', r2 );
	fprintf( 'x2=\n' );
	disp(x2);
	fprintf( 'r3=%e\n', r3 );
	fprintf( 'x3=\n' );
	disp(x3);

	fprintf( 'r%d=%e\n', n, r3 );

	% restart for better numerical condition
	if 0 == mod(n,10) && 0
		b = b - A*x3;
		x_appr = x_appr + x3;
	else
		x_appr = x3;
	end
end

%fprintf( 'x=\n' );
%disp(x);

fprintf( '|x-x_appr|=%e\n', norm(x-x_appr) );

