#ifndef ComplexLogF
#define ComplexLogF

#include <stdio.h>
#include <complex.h>


#define CMPLXE(x,y) ((x)+I*(y))

typedef struct COMPLEX {
	complex c;
} TypeComplex;

typedef struct COMPLEX_KAHAN {
	TypeComplex sum, c, cc, cs, ccs;
} TypeComplexKahan;

void initComplexKahan(TypeComplexKahan *ka);
void sumComplexKahan(TypeComplex x, TypeComplexKahan *ka);
TypeComplex totalComplexKahan(TypeComplexKahan *ka);

TypeComplex getComplex(double re, double im);
complex toLogComplex(TypeComplex a);
TypeComplex getNullComplex(void);
TypeComplex addComplex(TypeComplex a, TypeComplex b);
TypeComplex prodComplex(TypeComplex a, TypeComplex b);
TypeComplex divComplex(TypeComplex a, TypeComplex b);
TypeComplex conjComplex(TypeComplex a);
TypeComplex oppComplex(TypeComplex a);
double modComplex(TypeComplex a);
void fprintComplex(FILE *f, TypeComplex a);

double getRealComplex(TypeComplex a);
double getImagComplex(TypeComplex a);
double getRealLogComplex(TypeComplex a);
double getImagLogComplex(TypeComplex a);
#endif