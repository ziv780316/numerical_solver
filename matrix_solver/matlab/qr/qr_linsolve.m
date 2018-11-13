
function x = qr_linsolve( A, b )

	size_a = size(A);
	na = size_a(2);

	[Q R]=qr_householder(A', 1, 0);

	size_r = size(R);
	nr = size_r(2);

	if nr < na
		b(nr+1:end) = [];
	end

	C = R' * R;
	y = C\b;
	x = Q * R * y;

end
