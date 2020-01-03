#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 6;

double x0[] = {
	0.0, 
	0.0,
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
	double vdi = x[4];
	double vs = 0;
	double vsi = x[5];
	double vdisi = vdi - vsi;
	double vgs = vg - vs;
	double ids = (vgs - vth)*(vgs - vth) * (1 + vdisi*lamda);
	double gms = 2 * (vgs - vth) * (1 + vdisi*lamda);
	double gds = (vgs - vth) * (vgs - vth) * lamda;
	f[0] = vd - vdd; // KVL branch id
	f[1] = vg - vgg; // KVL branch ig
	f[2] = ig; // g
	f[3] = (vd - vdi) * gd + id; // d
	f[4] = (vdi - vd) * gd + ids; // di
	f[5] = (vsi - vs) * gs - ids; // si

	printf( "============= Circuit Information ============\n" );
	printf( "ids=%.15le\n", ids );
	printf( "gms=%.15le\n", gms );
	printf( "gds=%.15le\n", gds );
	printf( "vdisi=%.15le\n", vdisi );
}


void load_jacobian ( double *x, double *J )
{
}

bool bypass_check ( double *x, double *f, double *dx ) 
{
	bool bypass_violate = true;
	return bypass_violate;
}

