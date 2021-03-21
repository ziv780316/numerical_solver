function [c, A] = dft( f, T, m, n )
% ------------------------------------------
% dft perform Discrete-Fourier-Transform
% * INPUT:
% f = fitting function
% T = period
% m = Fourier harmonic order
% n = N sampling
%
% * OUTPUT:
% c = [c₋ₖ ... c₁ c₀ c₁ ... cₖ]ᵀ
% cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] cause ∂E(cₖ)/∂cₖ = 0
% A is c = A * f(t)
% ------------------------------------------
%
% * Mathematics background:
% Fourier curve fiting to minimize
% E(c) = 0.5 * ∑ⁿᵢ₌₁[f(tᵢ)-Pₘ(tᵢ)]
% Pₘ(t) = ∑ᵐₖ₌₋ₖ[cₖ*exp(j*2π*k*t)] -> two side
% Pₘ(t) = ∑²ᵐₖ₌₀[cₖ*exp(j*2π*k*t)] -> two side with mirrow central at m*ω (only valid when (2m+1) = n)
% DFT cₖ coefficent results is ∂E(cₖ)/∂cₖ = 0 for k = -m to m
% cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] cause ∂E(cₖ)/∂cₖ = 0
% MATLAB use cₖ = cₖ * n for normalize
%
% ------------------------------------------
%
% * DFT procedure:
% cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] is numerical integration (trapezoidal) of 
% cₖ = (1/T) * ∫ᵀ₀[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)dt]
% if f(t) ∈ ℜ then 
% cₖ = conj(c₋ₖ)
% computation complextity -> O(m*n) for c = A * f(t)
%
% ------------------------------------------
%
% * Practice consideration:
% if f(t) ∈ ℜ 
% if (2*m+1) = n then A⁻¹ exist and then has f(t) = (A⁻¹ * c) (∵ A is orthogonal matrix ∴ A⁻¹ = A*)
% [c₀ c₁ ... cₖ]ᵀ = A * [f(t₁) f(t₂) ... f(tₙ)]ᵀ
% Pₘ(t) = ∑ᵐₖ₌₋ₖ[cₖ*exp(j*2π*k*t)] has good to match f(t)
%
% if 2*(m + 1) > n
% then A⁻¹ does not exist ∴ Pₘ(t) = ∑ᵐₖ₌₋ₖ[cₖ*exp(j*2π*k*t)] is overdetermine to match f(t)
% there is 'alias' spectrum in the skirt ∴ Pₘ(t) is not accurate

% if 2*(m + 1) < n
% then A⁻¹ does not exist ∴ Pₘ(t) = ∑ᵐₖ₌₋ₖ[cₖ*exp(j*2π*k*t)] is underdetermine to match f(t)
%
% ------------------------------------------

if 4 == nargin
	real_f = false;
	debug = false;
elseif 5 == nargin
	debug = false;
end

t = linspace( 0, T, n+1 ); % 0, T/n, 2*(T/n), ..., T

x = transpose( f(t(1:n)) ); % f(0), f(T/n), ..., f((n-1)*(T/n))

A = zeros( 2*m+1, n ); 
for k = -m:1:m
	for i = 1:1:n
		A(k+m+1,i) = exp(-1j*2*pi*k*(i-1)/n);
	end
end

c = (1/n)*A*x; % [c₋ₖ ... c₁ c₀ c₁ ... cₖ]ᵀ

end % end of function dft
