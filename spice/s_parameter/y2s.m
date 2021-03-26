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
% E = [E0; 0] means feed in stimulus V = E0 in port 1
% i = E / (-Zo - I/Y) = (-Y*E) / (I + Y*Zo)
% v = -i / Y = E / (I + Y*Zo)
% V⁻ = (v + Zo * i)*0.5 
% V⁺ = (v - Zo * i)*0.5 = (E0/2) = Vav in stimulus port
% S = V⁻ / (E0/2) = (v + Zo * i) / E0
%
% ------------------------------------------
% S-parameter related to power:
%
% S₁₁ = (Zin₁ - Zo₁) / (Zin₁ + Zo₁) = V₁₁⁻ / Vav
% S₂₁ = V₂/(E/2) = V₂₂⁻ / Vav
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

E0 = 1;
n = size(Y, 1);
I = eye( n );
E = E0 * I;
if 1 == nargin 
	Zo = 50 .* I;
end
if nargin < 3
	debug_power = false;
end

i = (-Y*E) / (I + Y*Zo);
v = E / (I + Y*Zo);

%S = (I - Zo*Y) / (I + Zo*Y);
S = (v + Zo * i) / E0;

if debug_power && (n > 1)
	% analysis feed stimulus on port 1 
	S11 = S(1,1);
	S21 = S(2,1);
	Pav = ((E0*E0)/8)/Zo(1);
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
	V_neg = (v + Zo * i)*0.5;
	V_pos = (v - Zo * i)*0.5;
	fprintf( 'Vneg =\n' );
	disp( V_neg );
	fprintf( 'Vpos =\n' );
	disp( V_pos );
	P_Vneg = real(conj(V_neg) .* (V_neg/Zo) * 0.5);
	fprintf( 'P_Vneg =\n' );
	disp( P_Vneg );
	P_Vpos = real(conj(V_pos) .* (V_pos/Zo) * 0.5);
	fprintf( 'P_Vpos =\n' );
	disp( P_Vpos );
end

end % end of function y2s

