#include <math.h>
#include <limits.h>
#include "ComplexLog.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


static double complex clog1p(double complex a);
static double complex clogSumLog(double complex a, double complex b);


/* Kahan operations for complex numbers */

void initComplexKahan(TypeComplexKahan *ka) {
	ka->sum = getNullComplex();
	ka->c = getNullComplex();
	ka->cs = getNullComplex();
	ka->cc = getNullComplex();
	ka->ccs = getNullComplex();
}

void sumComplexKahan(TypeComplex x, TypeComplexKahan *ka) {
	TypeComplex t;
	t = addComplex(ka->sum, x);
	if(creal(ka->sum.c) >= creal(x.c))
		ka->c = addComplex(addComplex(ka->sum, oppComplex(t)), x);
	else
		ka->c = addComplex(addComplex(x, oppComplex(t)), ka->sum);
	ka->sum = t;
	t = addComplex(ka->cs, ka->c);
	if(creal(ka->cs.c) >= creal(ka->c.c))
		ka->cc = addComplex(addComplex(ka->cs, oppComplex(t)), ka->c);
	else
		ka->cc = addComplex(addComplex(ka->c, oppComplex(t)), ka->cs);
	ka->cs = t;
	ka->ccs = addComplex(ka->ccs, ka->cc);
}



TypeComplex totalComplexKahan(TypeComplexKahan *ka) {
	return addComplex(addComplex(ka->sum, ka->cs), ka->ccs);
}

/* Basic complex operations */

TypeComplex getComplex(double re, double im) {
  TypeComplex c;
  c.c = clog(CMPLXE(re, im));
  return c;
}

double complex toLogComplex(TypeComplex a) {
	return CMPLXE(creal(a.c)+log(fabs(cos(cimag(a.c)))), creal(a.c)+log(fabs(sin(cimag(a.c)))));
}

TypeComplex getNullComplex(void) {
  TypeComplex c;
  c.c = clog(CMPLXE(0., 0.));
  return c;
}

double getRealComplex(TypeComplex a) {
  return creal(cexp(a.c));
}

double getImagComplex(TypeComplex a) {
  return cimag(cexp(a.c));
}

double getRealLogComplex(TypeComplex a) {
  return creal((a.c));
}

double getImagLogComplex(TypeComplex a) {
  return cimag((a.c));
}

TypeComplex addComplex(TypeComplex a, TypeComplex b) {
  TypeComplex c;
  c.c = clogSumLog(a.c, b.c);
  return c;
}

TypeComplex prodComplex(TypeComplex a, TypeComplex b) {
  TypeComplex c;
  c.c = a.c+b.c;
  return c;
}

TypeComplex divComplex(TypeComplex a, TypeComplex b) {
  TypeComplex c;
  c.c = a.c-b.c;
  return c;
}

TypeComplex conjComplex(TypeComplex a) {
  TypeComplex c;
  c.c = CMPLXE(creal(a.c), -cimag(a.c));
  return c;
}

TypeComplex oppComplex(TypeComplex a) {
  TypeComplex c;
  c.c = CMPLXE(creal(a.c), cimag(a.c)+M_PI);
 return c;
}

double modComplex(TypeComplex a) {
  return exp(creal(a.c));
}

void fprintComplex(FILE *f, TypeComplex a) {
	TypeComplex c;
	c.c = cexp(a.c);
	if(cimag(c.c)>0.)
		fprintf(f, "%le\t+i\t%le\n", creal(c.c), cimag(c.c));
	if(cimag(c.c)<0.)
		fprintf(f, "%le\t-i\t%le\n", creal(c.c), fabs(cimag(c.c)));
	if(cimag(c.c)==0.)
		fprintf(f, "%le\t\n", creal(c.c));
}

double complex clog1p(double complex a) {
	if(cabs(a) < 1.E-07)
		return a;
	else
		return clog(1.0+a);
}

double complex clogSumLog(double complex a, double complex b) {
	double complex max, min;
	if(isinf(creal(a)) && creal(a) < 0)
		return b;
	if(isinf(creal(b)) && creal(b) < 0)
		return a;
	if(creal(a)>creal(b)) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	return max + clog1p(cexp(min-max));
}
