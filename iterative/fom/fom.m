%
rng(2);

use_rand = 1;
m = 4;

A = [...
1 -2 -2 1;...
10  3  1 5;...
5  1  2 -1;...
7  -3  3 -7;...
]; 

b = [...
4;...
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
maxiter = 10;

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

	if i > maxiter
		return;
	end

	% krylov subspace order
	n = i;

	% precondition
	A = P \ A;
	b = P \ b;

	% A*Q1 = Q2*H
	[Q1, Q2, H] = arnoldi_mgs( A, b, n );

	% fom
	Hn = Q1'*A*Q1;
	y = Hn \ (Q1'*b);
	y2 = (Hn \ eye(n,1)) * norm(b,2);
	z = Hn * y;
	x = Q1 * y;
	q_next = Q2(:,n+1);
	h_next = H(n+1,n);
	err = h_next * y(n,1) * q_next;
	b_proj = Q1*Hn*y;
	r = b - A*x;
	r2 = b - b_proj - err;
	r3 = -h_next * y(n,1) * q_next;
	%r
	%r2
	%r3
	%q_next
	fprintf( 'hnext=%.15e\n', h_next );
	fprintf( 'yn=%.15e\n', y(n,1) );
	fprintf( '|b - b_proj|=%.15e\n', norm(b - b_proj) );
	fprintf( '|r|=%.15e\n', norm(r) );
	%fprintf( '|r2|=%.15e\n', norm(r2) );
	%fprintf( '|err|=%.15e\n', norm(err) );
	fprintf( '|Q1 dot r|=%.15e\n', norm(Q1'*r,inf) );
end
