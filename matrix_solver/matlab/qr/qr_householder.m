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

function [Q, R, P, rank] = qr_householder( A, remove_redundant, debug )

	mat_size = size(A);
	nr = mat_size(1);
	nc = mat_size(2);
	
	min_norm = 1e-10;
	P = eye(nc);
	Ai = A;
	for i = 1:1:(nc - 1)
		v_norm = 0;
		pivot = i;
		while v_norm < min_norm
			ai = Ai(:,pivot);
			vi = ai;
			if i > 1
				vi(1:i-1) = 0;
			end

			v_norm = norm(vi, 2);
			% Ai is linear depend on other vector is v_norm -> 0
			if v_norm < min_norm
				pivot = pivot + 1;
			end
			if pivot > nc
				break;
			end
		end

		if pivot ~= i
			tmp = P(:,i);
			P(:,i) = P(:,pivot);
			P(:,pivot) = tmp;

			tmp = Ai(:,i);
			Ai(:,i) = Ai(:,pivot);
			Ai(:,pivot) = tmp;
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
			fprintf('Ai:\n');
			disp(Ai);
			fprintf('vi:\n');
			disp(vi);
			fprintf('wi:\n');
			disp(wi);
			fprintf('ui:\n');
			disp(ui);
			fprintf('Hi:\n');
			disp(Hi);
		end

		Ai = Hi * Ai;
	end

	R = Ai;

	% count rank
	rank = 0;
	for i = 1:1:nc
		if abs(R(i,i)) > 1e-10
			rank = rank + 1;
		end
	end

	if remove_redundant
		k = 1;
		for i = 1:1:nc
			if abs(R(k,i)) < 1e-10
				fprintf('delete redundant column %d\n', i);
				R(k,:) = [];
			end
			k = k + 1;
		end
	end

end
