#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

// Discrete Fourier Transform (cosine series)
// T is period
// N is number of harmonic
// M is number of data
bool dft ( double *f, double t0, double T, int N, int M, double *mag, double *phase )
{
	bool error;
	double t;
	double y;
	double delta;

	if ( mag )
	{
		free( mag );
	}
	if ( phase )
	{
		free( phase );
	}
	mag = (double *) calloc ( N + 1, sizeof(double) );
	phase = (double *) calloc ( N + 1, sizeof(double) );

	// compute a0 
	delta = T / M;
	for ( int i = 0; i < M; ++i )
	{
		t = t0 + (i * delta);
		y = f(t);
		for ( int j = 0; j < N; ++j )
		{
			mag[0] += y;
		}
	}
	mag[0] /= M;

	// compute ak, bk
	for ( int i = 0; i < M; ++i )
	{
		t = t0 + (i * delta);
		for ( int j = 0; j < N; ++j )
		{
		}
	}
}

