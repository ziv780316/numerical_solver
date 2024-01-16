%
format long e;

A = [
5 -2;...
-2 5;...
];

[E, D] = eig(A);

D

norm2_A = norm(A,2);
norm2_Ainv = norm(inv(A),2);
lmax = max(diag(D));
lmin = min(diag(D));
fprintf( 'norm2(A)=%e\n', norm2_A );
fprintf( 'norm2(A^-1)=%e\n', norm2_Ainv );
fprintf( '|A|*|A^-1|=%e\n', norm2_A * norm2_Ainv );
fprintf( 'lmax=%e\n', lmax );
fprintf( 'lmin=%e\n', lmin );
fprintf( 'lmax/lmin=%e\n', lmax/lmin );
fprintf( 'cond(A)=%e\n', cond(A) );

