#ifndef HOMOTOPY_H
#define HOMOTOPY_H

typedef enum {
	HOMOTOPY_EXTRAPOLATE_NONE,
	HOMOTOPY_EXTRAPOLATE_DIFFERENCE,
	HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL,
} homotopy_extrapolate_type;

typedef struct
{
	homotopy_extrapolate_type extrapolate_type;

	double lamda_start;
	double lamda_stop;

} homotopy_param_t;


#endif

