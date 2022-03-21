rng(1);
m = 3;
A = [...
0 -2 -2;...
1  3  1;...
0  0  2;...
]; % minimum polynomial degree is 2
%m = minpoly(A);
b = rand(m,1);

D = diag(diag(A));
P = eye(m,m);
%P = D;
%P = full(ilu(sparse(A)));

A = P\A;
b = P\b;
n = 2;
[Q1, Q2, H] = arnoldi_mgs( A, b, n );

H_sub = H(1:n,:);
x = A\b;
y = H_sub\(Q1'*b);
x_appr = Q1*y;

fprintf( 'x=\n' );
disp(x);
fprintf( 'x_appr=\n' );
disp(x_appr);

