%
format long e;

% Fibonacci LFSR x^4 + x^3 + 1
A = [...
0 1 0 0;...
0 0 1 0;...
0 0 0 1;...
1 1 0 0;... % taps bit x^4 + x^3 + 1
];

C = transpose(A);

% s0(1) is LSB, s0(n) is MSB
s0 = [...
1;...
1;...
1;...
1;...
];

% period = 2^n - 1
period = 15;
for i = 1:1:period
	sk = get_s_next( A, s0 );
	fprintf( 'sk_%d=\n', i );
	disp( sk );
	s0 = sk;
end

% A^period willl be I
Ak = get_Ak( A, period );
fprintf( 'Ak_%d=\n', i );
disp( Ak );

function sk = get_s_next( A, s0 )
	sk = A * s0;
	n = size(s0);
	for i = 1:1:n
		if sk(i) > 1
			sk(i) = 0;
		end
	end
end 

function Ak = get_Ak( A, k )
	n = size(A);
	for m = 1:1:k
		if m == 1
			Ak = A;
		else
			Ak = Ak * A;
		end

		for i = 1:1:n
			for j = 1:1:n
				if Ak(i,j) > 1
					Ak(i,j) = 0;
				end
			end
		end
	end
end 
