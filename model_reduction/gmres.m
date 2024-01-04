%
rng(1);

use_rand = 1;
m = 4;

A = [...
1 -2 -2 1;...
1  3  1 5;...
5  1  2 -1;...
7  -3  3 -7;...
]; 

b = [...
1;...
2;...
3;...
4;...
];

if use_rand
	m = 10;
	A = rand(m,m);
	b = rand(m,1);
end

fprintf( 'rank(A)=%d\n', rank(A) );

debug = 0;

%preconditioner = 'no';
%preconditioner = 'diagnoal';
preconditioner = 'ilu';
if strcmp(preconditioner,'no')
	P = eye(m,m);
elseif strcmp(preconditioner,'diagonal')
	P = diag(diag(A));
elseif strcmp(preconditioner,'ilu')
	P = full(ilu(sparse(A)));
end

for i = 1:1:m
	fprintf( 'iter=%d\n', i);

	% krylov subspace order
	n = i;

	% precondition
	A = P \ A;
	b = P \ b;

	% A*Q1 = Q2*H
	[Q1, Q2, H] = arnoldi_mgs( A, b, n );

	% gmres
	AQ = A*Q1;
	%y = (AQ'*AQ)\(AQ'*b);
	y = (H'*H)\(H'*(norm(b))*eye(n+1,1));
	x = Q1 * y;
	r = b - A*x;
	fprintf( '|r|=%.15e\n', norm(r) );
	fprintf( '|AQ dot r|=%.15e\n', norm(AQ'*r,inf) );

	b_proj = AQ*y;
	b_proj_proj_AQ = AQ*inv((AQ'*AQ))*AQ'*b_proj;

	if debug
		fprintf( 'b_proj=\n' );
		disp( b_proj );
		fprintf( 'b_proj_proj_AQ=\n' );
		disp( b_proj_proj_AQ );
	end

	fprintf( '|b_proj - b_proj_proj_AQ|=%.15e\n', norm(b_proj-b_proj_proj_AQ, inf) );

	if debug
		fprintf( 'b_proj=\n' );
		disp(b_proj);
		fprintf( 'r=\n' );
		disp(r);
		fprintf( 'Q1*r=\n' );
		disp(Q1'*r);
	end
end
