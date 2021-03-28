% lossless network
clear;
format long e;
f = 1/(2*pi);
c = 1;
Y = [1j*2*pi*f*c  -1j*2*pi*f*c; -1j*2*pi*f*c 1j*2*pi*f*c];
Zo = 1 * eye(2);
[S,v,i] = y2s(Y, Zo, true);
fprintf( 'S=\n' );
disp( S );
fprintf( 'v=\n' );
disp( v );
fprintf( 'i=\n' );
disp( i );

