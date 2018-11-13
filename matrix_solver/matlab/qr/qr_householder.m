% QR decomposition by Householder reflection
% A  = [a1 a2 ... an]
% Ai = [a1i a2i ... ani] --> A of ith iteration
% vi = [aii] and vi(1:i-1) = 0 if i - 1 > 0
% wi = [0 0 ... 0] with wi(i) = |vi| --> norm2
% ui = (vi - wi) / norm2(vi - wi) 
% Hi = I - 2 * u * u' --> Householder reflection
% Ai+1 = Hi * Ai
% Hn * ... * H1 * A = R --> n = col - 1
% Q = H1' * ... * Hn'

function [Q, R, rank] = qr_householder( A, remove_redundant, debug )

	mat_size = size(A);
	nr = mat_size(1);
	nc = mat_size(2);
	rank = nc;
	
	Ai = A;
	for i = 1:1:(nc - 1)
		vi = Ai(:, i);
		if i > 1
			vi(1:i-1) = 0;
		end
		wi = zeros(nr, 1);
		wi(i) = norm(vi, 2);
		ui = (vi - wi) / norm(vi - wi, 2);
		Hi = eye(nr) - 2*ui*ui';
		if 1 == i
			Q = Hi';
		else
			Q = Q * Hi';
		end

		if debug
			disp(Ai);
			disp(vi);
			disp(wi);
			disp(ui);
			disp(Hi);
		end

		Ai = Hi * Ai;
	end

	R = Ai;

	if remove_redundant
		for i = 1:1:nc
			if abs(R(i,i)) < 1e-12
				fprintf('delete redundant column %d\n', i);
				R(:,i) = [];
			end
		end
	end

end
