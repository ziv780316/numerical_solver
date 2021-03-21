% reset
clear;
close all;
format long e;

f1 = @(t) 1 + sin(2*pi*1*t) + 3*sin(2*pi*3*t);
T = 1;
m = 9;
n = 19;
m_matlab = floor(n/2);
t = linspace(0, T, n+1);
[c, A] = dft( f1, T, m, n );
X = fft( f1(t), n );
c_matlab = transpose((1/n) * fftshift( X ));
for k = -m:1:m
	fprintf( 'c%d=%.10e + j*%.10e, |c%d|=%.10e \n', k, real(c(k+m+1)), imag(c(k+m+1)), k, abs(c(k+m+1)) );
end
for k = -m_matlab:1:m_matlab
	fprintf( 'c_matlab_%d=%.10e + j*%.10e, |c%d|=%.10e \n', k, real(c_matlab(k+m_matlab+1)), imag(c_matlab(k+m_matlab+1)), k, abs(c_matlab(k+m_matlab+1)) );
end

subplot(2, 2, 1);
hold on;
t = linspace(0, T, 100*n);
plot( t, f1(t), 'b-' );

subplot(2, 2, 2);
hold on;
plot( t, f1(t), 'b-' );

t = linspace(0, T, n+1);
[x, A_inv] = idft(c, T, m, t );
subplot(2, 2, 1);
plot( t, x', 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r' );

t = linspace(0, T, n+1);
[x, A_inv] = idft(c_matlab, T, m_matlab, t );
subplot(2, 2, 2);
plot( t, x', 'g-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'g' );

subplot(2, 2, 3);
dft_spectrum( c, T );

subplot(2, 2, 4);
dft_spectrum( c_matlab, T );

