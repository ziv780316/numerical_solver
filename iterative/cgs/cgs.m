%
rng(2);

use_rand = 1;

%m = 4;
%
%A = [...
%1 -2 -2 1;...
%10  3  1 5;...
%5  1  2 -1;...
%7  -3  3 -7;...
%]; 

%b = [...
%4;...
%2;...
%3;...
%4;...
%];

use_rand = 0;

m = 3;

A = [...
2 1 1;...
1 2 1;...
1 1 2;...
]; 

b = [...
4;...
2;...
3;...
];

if use_rand
	m = 10;
	A = rand(m,m);
	b = rand(m,1);
	A = A'*A;
	%A = A + 2*eye(m,m);
end


fprintf( 'rank(A)=%d\n', rank(A) );
fprintf( 'cond(A)=%d\n', cond(A) );

debug = 0;
maxiter = 3;

x = zeros(m,1);
r = b - A*x;
r0 = r;
e = 0.5*x'*A*x - x'*b;
p = r;
p0 = p;
P = zeros(m,m);
R = zeros(m,m);
a = zeros(m,1);
fprintf( '|r0|=%.15e\n', norm(r) );
fprintf( '|e0|=%.15e\n', e );
cnt = 0;

for i = 1:1:m
	fprintf( 'iter=%d\n', i);

	if i > maxiter
		break;
	end


	alpha = (r'*r) / (p'*A*p);
	a(i) = alpha;
	P(:,i) = p;
	R(:,i) = r;
	cnt = cnt + 1;

	% krylov subspace order
	n = i;
	[Q, Q2, H] = arnoldi_mgs( A, b, n );
	b_proj_Q = Q*inv((Q'*Q))*Q'*b;
	Pi=P(:,1:n);
	b_proj_P = Pi*inv((Pi'*Pi))*Pi'*b;
	fprintf( 'b_proj_Q=\n' );
	disp( b_proj_Q );
	fprintf( 'b_proj_P=\n' );
	disp( b_proj_P );
	fprintf( 'b_proj_P - b_proj_Q=\n' );
	disp( b_proj_P - b_proj_Q );


	x = x + alpha * p;
	e = 0.5*x'*A*x - x'*b;
	r_old = r;
	r = r_old - alpha * A * p;
	%r = b - A*x;
	
	beta = (r'*r) / (r_old'*r_old);
	p_old = p;
	p = r + beta * p;


	fprintf( 'alpha=%.15e\n', alpha );
	fprintf( 'beta=%.15e\n', beta );
	fprintf( '|r|=%.15e\n', norm(r) );
	fprintf( '|e|=%.15e\n', e );
	fprintf( '<r,p>=%.15e\n', r'*p );
	fprintf( '<r,p0>=%.15e\n', r'*p0 );
	fprintf( '<r,p_old>=%.15e\n', r'*p_old );
	fprintf( '<r,r0>=%.15e\n', r'*r0 );
	fprintf( '<r,r_old>=%.15e\n', r'*r_old );
end
x_check = P(:,1:cnt) * a(1:cnt,1);
%x
%x_check
%x - x_check

P'*A*P*a - P'*b
