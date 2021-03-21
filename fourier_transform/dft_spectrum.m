function dft_spectrum( c, fs );
% Plot DFT spectrum
% * INPUT:
% c = cₖ, k = -m to m, cₖ = (1/n) * ∑ⁿᵢ₌₁[f(tᵢ)*exp(-j*2π*k*ω*tᵢ)] 
% T = period
% ------------------------------------------

m = floor(max(size(c))/2);
k = -m:1:m;
f = k * fs;

figure;
stem( f, abs(c)', 'bo', 'MarkerSize', 6 );
grid on;
xlabel( 'Frequency (Hz)' );
ylabel( 'Amplitude (Mag)' );
title( 'Spectrum' );
end % end of function dft_spectrum
