function [S,v,i] = y2s( Y, Zo, debug_power )
% y2s
% ------------------------------------------
% Mathematics background:
%  
% /-------\ /---\   /---\
% | Y | I | | v |   | 0 |  
% |---|---------- = -----  
% | I |-Zo| | i |   | E |
% \-------/ \---/   \---/
%
% E = [1;0] means feed in stimulus V=1 in port 1
% i = E / (-Zo - I/Y) = (-Y*E) / (I + Y*Zo)
% v = -i / Y = E / (I + Y*Zo)
% S = v + Zo * i
%
% ------------------------------------------
% S-parameter related to power:
%
% S₁₁ = (Zin₁ - Zo₁) / (Zin₁ + Zo₁)
% S₂₁ = V₂/(E/2)
% |S₁₁|²= (Pav - P₁)/Pav
% |S₂₁|²= (P₂)/Pav
% 
% lossless network P₁ = P₂ thus |S₁₁|² + |S₂₁|² = (Pav + P₂ - P₁)/Pav = 1
% reciprocal network S₁₁ = S₂₂ and S₂₁ = S₁₂
% 
% Sᴴ*S = 
% /----------------------------------------\
% |    |S₁₁|²+|S₂₁|²   | S₁₁'S₁₂ + S₂₁'S₂₂ |
% |--------------------|-------------------| = I
% | S₁₁'S₁₂ + S₂₂'S₂₁  |   |S₂₂|²+|S₁₂|²   |
% \----------------------------------------/
% ------------------------------------------

n = size(Y, 1);
I = eye( n );
if 1 == nargin 
	Zo = 50 .* I;
end
if nargin < 3
	debug_power = false;
end

i = (-Y*I) / (I + Y*Zo);
v = I / (I + Y*Zo);

S = (I - Zo*Y) / (I + Zo*Y);

if debug_power && (n > 1)
	% analysis feed stimulus on port 1 
	S11 = S(1,1);
	S21 = S(2,1);
	Pav = (1/8)/Zo(1);
	is = -i(1,1);
	P1 = real(conj(v(1,1))*is*0.5);
	P2 = real(conj(v(2,1))*is*0.5);
	P1_s = abs((conj(S11)*S11)*Pav-Pav);
	P2_s = abs((conj(S21)*S21)*Pav);
	Ploss = abs(((conj(S11)*S11 + conj(S21)*S21)*Pav - Pav));
	vz1 = is*Zo(1);
	Pz1 = real(conj(vz1)*is*0.5);
	fprintf( 'Pav = %.10e\n', Pav );
	fprintf( 'Preturn = %.10e\n', Pav-P1 );
	fprintf( 'Pz1 = %.10e\n', Pz1 );
	fprintf( 'P1 = %.10e\n', P1 );
	fprintf( 'P2 = %.10e\n', P2 );
	fprintf( 'P1_s = %.10e\n', P1_s );
	fprintf( 'P2_s = %.10e\n', P2_s );
	fprintf( 'P1-P2 = %.10e\n', P1-P2 );
	fprintf( 'Ploss = %.10e\n', Ploss );
end

end % end of function y2s

