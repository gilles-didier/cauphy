#include <limits.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "SimulateTrait.h"

static void simulTraitRec(int n, double anc, double *trait, TypeTree *tree, TypeFunctDrawTrans trans, void *data);

double cauchyInitStd(void *data) {
	return ((TypeCauchyStdSimulData*)data)->disp*((TypeCauchyStdSimulData*)data)->bl*tan(M_PI*(unif_rand()-0.5));
}

double cauchyInitStart(void *data) {
	return ((TypeCauchyStartSimulData*)data)->start;
}

double cauchyOUInit(void *data) {
	return ((TypeOUCauchyParam*)data)->start;
}

double cauchyOUTrans(double val, double time, void *data) {
	return (((TypeOUCauchyParam*)data)->disp*(1-exp(-((TypeOUCauchyParam*)data)->str*time))/((TypeOUCauchyParam*)data)->str)*tan(M_PI*(unif_rand()-0.5))+exp(-((TypeOUCauchyParam*)data)->str*time)*val+((TypeOUCauchyParam*)data)->mean*(1.-exp(-((TypeOUCauchyParam*)data)->str*time));
}

double cauchyTransStd(double val, double time, void *data) {
	return val+((TypeCauchyStdSimulData*)data)->disp*time*tan(M_PI*(unif_rand()-0.5));
}

double cauchyTransStart(double val, double time, void *data) {
	return val+((TypeCauchyStartSimulData*)data)->disp*time*tan(M_PI*(unif_rand()-0.5));
}

void simulTraitRec(int n, double anc, double *trait, TypeTree *tree, TypeFunctDrawTrans trans, void *data) {
	int c;
	trait[n] = trans(anc, tree->time[n], data);
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		simulTraitRec(c, trait[n], trait, tree, trans, data);
}

void simulTrait(double *trait, TypeTree *tree, TypeFunctDrawInit init, TypeFunctDrawTrans trans, void *data) {
	int c;
	trait[tree->root] = init(data);
	for(c=tree->node[tree->root].child; c!=NOSUCH; c=tree->node[c].sibling)
		simulTraitRec(c, trait[tree->root], trait, tree, trans, data);
}
