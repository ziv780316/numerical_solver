% Compute eigenvalue by QR decomposition 
% Ai = Qi*Ri
% Ri = Qi^(-1)*A 
% Ai+i = Ri*Qi = Qi^(-1)*Ai*Q  --> Ai+1 and Ai is similar matrix (the same eigenvalue)
% when i -> inf, diag(Ri) is eigenvalues of A

function [E,V] = qr_eig( A, debug )

	tol = 1e-10;
	size_A = size(A);
	nrow = size_A(1);
	ncol = size_A(2);

	E = inf .* eye(ncol);
	iter = 1;
	while true
	
		[Q,R,P,rank] = qr_householder(A,0,0);

		if debug
			fprintf('i=%d\n',iter);
			fprintf('A:\n');
			disp(A);
			fprintf('Q:\n');
			disp(Q);
			fprintf('R:\n');
			disp(R);
		end

		A = R*Q;
		converge = true;
		for i = 1:1:ncol
			if abs(A(i,i) - E(i,i)) > tol
				converge = false;
				break;
			end
		end

		E = A;
		if converge
			break;
		else
			iter = iter + 1;
		end
	end

	E = diag( A );

	V = zeros(ncol,ncol);
	
end
