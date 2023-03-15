#ifndef CauchyF
#define CauchyF
#include "Tree.h"
#include "ComplexLog.h"

typedef struct CAUCHY_PARAM {
	double start, disp;
} TypeCauchyParam;

typedef struct CAUCHY_INFO {
	int sizeChild, *child;
	double *time;
	TypeComplex *A;
} TypeCauchyInfo;


#ifdef __cplusplus
extern "C" {
namespace Tree {
#endif

double getCauchyLogDensityStandard(double x,  double d);
void fillCauchyInfo(int n, TypeTree *tree, double param, TypeCauchyInfo *cinf);
void freeCauchyInfo(int n, TypeTree *tree, TypeCauchyInfo *cinf);
double getCauchyLogDensityStem(TypeCauchyInfo cinf, double *val, double start);
double getCauchyLogDensityNoStem(TypeCauchyInfo cinfL, TypeCauchyInfo cinfR, double *val, double start);
void fillCauchyAncestralPosteriorLogDensityNoStem(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start);
void fillCauchyAncestralPosteriorLogDensityStem(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start);
void fillCauchyAncestralPosteriorLogDensityREML(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp);
void fillCauchyIncrementPosteriorLogDensityNoStem(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start);
void fillCauchyIncrementPosteriorLogDensityStem(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double disp, double start);
void fillCauchyIncrementPosteriorLogDensityREML(int n, double *dens, double *tabVal, int nVal, TypeTree *tree, double dis);

#ifdef __cplusplus
}
}
#endif

#endif
