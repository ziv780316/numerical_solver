#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 4;

double x0[] = {
	0.0, 
	0.0,
	0.0,
	0.0,
};

// vd d 0 1
// vg g 0 0.7
// m1 d g 0 0 model_1
void load_f ( double *x, double *f )
{	
	double gd = 1;
	double gs = 1;
	double vth = 0.3;
	double lamda = 1e-3;
	double vgg = 0.7;
	double vdd = 1;
	double id = x[0];
	double ig = x[1];
	double vg = x[2];
	double vd = x[3];
	double vdi = vd - (-id) * gd;
	double vs = 0;
	double vsi = vs + -id * gs;
	double vdisi = vdi - vsi;
	double vgs = vg - vs;
	double ids = (vgs - vth)*(vgs - vth) * (1 + vdisi*lamda);
	double gms = 2 * (vgs - vth) * (1 + vdisi*lamda);
	double gds = (vgs - vth) * (vgs - vth) * lamda;
	double y11 = 0;
	double y12 = 0;
	double y1b = ids;
	double fdi;
	double fsi;
	f[0] = vd - vdd; // KVL branch id
	f[1] = vg - vgg; // KVL branch ig
	f[2] = ig; // g
	f[3] = id + (vd * y11) + (vs * y12) + y1b ; // d
	fdi = (vdi - vd) * gd + ids; // di
	fsi = (vsi - vs) * gs - ids; // si

	printf( "============= Circuit Information ============\n" );
	printf( "ids=%.15le\n", ids );
	printf( "gms=%.15le\n", gms );
	printf( "gds=%.15le\n", gds );
	printf( "vdi=%.15le\n", vdi );
	printf( "vsi=%.15le\n", vsi );
	printf( "vdisi=%.15le\n", vdisi );
	printf( "fdi=%.15le\n", fdi );
	printf( "fsi=%.15le\n", fsi );
}


void load_jacobian ( double *x, double *J )
{
}

bool bypass_check ( double *x, double *f, double *dx ) 
{
	bool bypass_violate = true;
	return bypass_violate;
}

