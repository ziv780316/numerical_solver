
% underdetermined circuit system
% v1 1 0 1
% v2 2 0 3
% v3 2 1 2

A = [ 0 0 1 -1 0;
      0 0 0 1 1;
      1 0 0 0 0;
      0 1 0 0 0;
      -1 1 0 0 0 ];

b = [0;0;1;3;2];

x = qr_linsolve(A,b,0);
r = A*x-b;

fprintf('x:\n');
disp(x);
fprintf('residual:\n');
disp(r);
