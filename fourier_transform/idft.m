function [x, A_inv] = idft( c, T, m, t )
% ------------------------------------------
% Perform polynomial interpolation by DFT results
% * INPUT:
% c = cₖ, k = -m to m, cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] 
% T = period
% m = Fourier harmonic order
% n = N sampling
% t = time vector to interpolate f(t)
%
% * OUTPUT:
% x is interpolated vector [Pₘ(t₁) ... Pₘ(tₙ)]ᵀ
% A_inv is A⁻¹ -> x = A⁻¹ * c
% ------------------------------------------
%
% * Mathematics background:
% Pₘ(t) = ∑ᵐₖ₌₋ₖ[cₖ*exp(j*2π*k*t)]
% cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] cause ∂E(cₖ)/∂cₖ = 0
%
% ------------------------------------------
%
% if f(t) ∈ ℜ then 
% cₖ = conj(c₋ₖ)
% Pₘ(t) = c₀ + 2*Re{∑ᵐₖ₌₁[cₖ*exp(j*2π*k*t)]}
%
% ------------------------------------------

t_len = max(size(t));

A_inv = zeros( t_len, 2*m+1 ); 
for i = 1:1:t_len
	for k = -m:1:m
		A_inv(i, k+m+1) = exp(j*2*pi*k*t(i)/T);
	end
end

x = A_inv * c;
x = real(x);

end % end of function idft
