rng(1);
m = 10;
A = rand(m,m) + eye(m,m);
b = rand(m,1);

D = diag(diag(A));
P = eye(m,m);
%P = D;
%P = full(ilu(sparse(A)));

A = P\A;
b = P\b;
n = 8;
[Q1, Q2, H] = arnoldi_mgs( A, b, n );

H_appr = H(1:n,:);
x = A\b;

% chord newton
iter_max = 10;
x_appr = zeros(m,1);
for i =1:1:iter_max
  f = A*x_appr - b;
  y = H_appr\(Q1'*f);
  dx = -Q1*y;
  %y = H\(Q2'*f);
  %dx = -Q1*y;
  x_appr = x_appr + dx;
  fprintf( '|f%d|=%e\n', i, norm(f) );
  fprintf( '|dx%d|=%e\n', i, norm(dx) );
end

fprintf( 'x=\n' );
disp(x);
fprintf( 'x_appr=\n' );
disp(x_appr);

