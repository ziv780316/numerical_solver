% -----------------------------------
% [netlist]
% p1 in1 0 z0=z01
% p2 in2 0 z0=z02
% r1 in1 in2 0 r
% -----------------------------------
% [matrix]
% node1 = in1
% node2 = in2
% node3 = p1:p
% node4 = p2:p
% -----------------------------------

format long e;

debug = false;
if ~isempty( getenv('DEBUG_COMPUTE_CZ') )
	debug = true;
	fprintf( 'turn of debug mode of compute_cz\n' );
end

% ------------------------------------------------
% input setup
k = 1.380649e-23;
T = 290;
t_coef = 4*k*T;

z01 = 50.0;
z02 = 50.0;
r = 10.0;
g = 1/r;

n = 4;

J = [...
 g -g    1    0;...
-g  g    0    1;...
 1  0 -z01    0;...
 0  1    0 -z02;...
];

Z0 = [...
z01   0;...
0   z02;...
];

if debug
	fprintf( 'N=%d\n', n );
	fprintf( 'J=\n' );
	disp( J );
	fprintf( 'Z0=\n' );
	disp( Z0 );
end

% ------------------------------------------------
% compute xf_2
rhs = zeros(n, 1);
rhs(1) = 1;
xf1 = transpose(J) \ rhs;
if debug
	fprintf( '-------------------------------\n', n );
	fprintf( 'compute xf1:\n' );
	fprintf( 'rhs=\n' );
	disp( rhs );
	fprintf( 'xf1=\n' );
	disp( xf1 );
end

rhs = zeros(n, 1);
rhs(2) = 1;
xf2 = transpose(J) \ rhs;

xf_n = [xf1, xf2];
if debug
	fprintf( '\n' );
	fprintf( 'compute xf2:\n' );
	fprintf( 'rhs=\n' );
	disp( rhs );
	fprintf( 'xf2=\n' );
	disp( xf2 );
end

xf_2 = [xf1(1:2), xf2(1:2)];

if debug
	fprintf( '\n' );
	fprintf( 'xf_n:\n' );
	disp( xf_n );
	fprintf( 'xf_2:\n' );
	disp( xf_2 );
end

% ------------------------------------------------
% compute cz
cz11 = t_coef * g * (xf1(1) - xf1(2))^2;
cz12 = t_coef * g * (xf1(1) - xf1(2)) * conj(xf2(1) - xf2(2));
cz21 = conj(cz12);
cz22 = t_coef * g * (xf2(1) - xf2(2))^2;
cz = [cz11, cz12; cz21, cz22] ./ t_coef;
if debug
	fprintf( '-------------------------------\n', n );
	fprintf( 'compute cz:\n' );
	fprintf( 'cz=\n' );
	disp( cz );
end

% ------------------------------------------------
% compute cy
cy = inv(xf_2') * cz * inv(xf_2);
if debug
	fprintf( '-------------------------------\n', n );
	fprintf( 'compute cy:\n' );
	fprintf( 'cy=\n' );
	disp( cy );
end
