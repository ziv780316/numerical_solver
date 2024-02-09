%
format long e;

m = 2;

A = [...
2 1;...
1 2;...
]; 

b = [...
1;...
3;...
];

x0 = [...
2;...
-4;...
];

plot_range = 4;

fprintf( 'rank(A)=%d\n', rank(A) );
fprintf( 'cond(A)=%d\n', cond(A) );

debug = 1;
maxiter = 2;

% conjugate gradient 
x = x0;
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
xp = [x(1)];
yp = [x(2)];
tol = 1e-10;

for i = 1:1:m
	fprintf( 'iter=%d\n', i);

	if i > maxiter
		break;
	end

	if norm(r) < tol 
		break;
	end

	alpha = (r'*r) / (p'*A*p);
	a(i) = alpha;
	P(:,i) = p;
	R(:,i) = r;
	cnt = cnt + 1;

	x = x + alpha * p;
	xp = [xp, x(1)];
	yp = [yp, x(2)];
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
	if debug
		fprintf( '<r,p>=%.15e\n', r'*p );
		fprintf( '<r,p0>=%.15e\n', r'*p0 );
		fprintf( '<r,p_old>=%.15e\n', r'*p_old );
		fprintf( '<r,r0>=%.15e\n', r'*r0 );
		fprintf( '<r,r_old>=%.15e\n', r'*r_old );
	end
end

fp = fcontour(@(x,y) 0.5*[x y]*A*[x;y] - [x y]*b, range.*[-1 1 -1 1]);
fp.LevelList = range.*linspace(0,10,11);
hold on;
plt = plot(xp, yp, 'b-pentagram', 'linewidth', 2, 'markersize', 15);
plt.MarkerFaceColor = [0 0 1];

% steepest descent
x = x0;
r = b - A*x;
r0 = r;
e = 0.5*x'*A*x - x'*b;
R = zeros(m,m);
fprintf( '|r0|=%.15e\n', norm(r) );
fprintf( '|e0|=%.15e\n', e );
cnt = 0;
xp = [x(1)];
yp = [x(2)];

maxiter = 10;

for i = 1:1:maxiter
	fprintf( 'iter=%d\n', i);

	if i > maxiter
		break;
	end

	if norm(r) < tol 
		break;
	end

	alpha = (r'*r) / (r'*A*r);
	R(:,i) = r;
	cnt = cnt + 1;

	x = x + alpha * r;
	xp = [xp, x(1)];
	yp = [yp, x(2)];
	e = 0.5*x'*A*x - x'*b;
	r_old = r;
	r = b - A*x;
	
	fprintf( '|r|=%.15e\n', norm(r) );
	fprintf( '|e|=%.15e\n', e );
	if debug
		fprintf( '<r,r_old>=%.15e\n', r'*r_old );
	end
end

plt = plot(xp, yp, 'r-', 'linewidth', 2, 'markersize', 15);
plt.MarkerFaceColor = [0 0 1];
