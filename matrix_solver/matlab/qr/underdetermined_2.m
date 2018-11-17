
% underdetermined circuit system
% v1 1 0 1
% i1 2 1 2
% i2 2 0 3

A = [ 0 0 1;
      0 0 0;
      1 0 0 ];

b = [2;-3;1];

x = qr_linsolve(A,b,0);
r = A*x-b;

fprintf('x:\n');
disp(x);
fprintf('residual:\n');
disp(r);
