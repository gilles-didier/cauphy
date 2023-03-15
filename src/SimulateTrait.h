#ifndef SimulateTraitF
#define SimulateTraitF
#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Tree.h"
#include "Cauchy.h"

typedef double TypeFunctDrawInit(void*);
typedef double TypeFunctDrawTrans(double, double, void*);

typedef struct CAUCHY_DATA_START {
	double disp, start;
} TypeCauchyStartSimulData;

typedef struct CAUCHY_DATA_STD {
	double disp, bl;
} TypeCauchyStdSimulData;

double cauchyInitStart(void *data);
double cauchyTransStart(double val, double time, void *data);
void simulTrait(double *trait, TypeTree *tree, TypeFunctDrawInit init, TypeFunctDrawTrans trans, void *data);

#endif
