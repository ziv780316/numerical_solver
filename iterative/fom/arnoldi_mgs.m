function [Q1, Q2, H] = arnoldi_mgs( A, b, n )
% Arnoldi-MGS decomposition
% A*Q1 = Q2*H

m = size(A,1);
Q2 = zeros(m, n+1);
H = zeros(n+1, n);

q = b / norm(b);
Q2(:,1) = q;
for i = 1:1:n
  q = Q2(:,i);
  u = A*q;
  for k = 1:1:i
    q = Q2(:,k);
    h = u'*q;
    u = u - h*q;
    H(k,i) = h;
  end
  h = norm(u);
  H(i+1,i) = h;
  u = u / h;
  Q2(:,i+1) = u;
end
Q1 = Q2(:,1:n);

end
