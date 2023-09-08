#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "Kahan.h"


void initDoubleKahan(TypeDoubleKahan *ka) {
	ka->sum = 0.;
	ka->c = 0.;
	ka->cs = 0.;
	ka->cc = 0.;
	ka->ccs = 0.;
}

void sumDoubleKahan(double x, TypeDoubleKahan *ka) {
	double t;
	t = ka->sum+x;
	if(fabs(ka->sum) >= fabs(x))
		ka->c = ka->sum-t+x;
	else
		ka->c = x-t+ka->sum;
	ka->sum = t;
	t = ka->cs+ka->c;
	if(fabs(ka->cs) >= fabs(ka->c))
		ka->cc = ka->cs-t+ka->c;
	else
		ka->cc = ka->c-t+ka->cs;
	ka->cs = t;
	ka->ccs = ka->ccs+ka->cc;
}

double totalDoubleKahan(TypeDoubleKahan *ka) {
	return ka->sum+ka->cs+ka->ccs;
}

double sumSignedLogTableKahan(double a[], int sign[], int size) {
	int i, npos, max_sign, max_ind;
	double max_double, res;
	if(size <= 0)
		return log(0.);
	max_ind = 0;
	max_double = a[max_ind];
	max_sign = sign[max_ind];
	for (i = 1; i < size; i++) {
		if (a[i] > max_double) {
			max_ind = i;
			max_double = a[i];
			max_sign = sign[i];
		}
	}
	npos = 0;
	for (i = 0; i < size; i++)
		if (sign[i] >= 1)
			npos++;
	TypeDoubleKahan kexpSum;
	initDoubleKahan(&kexpSum);
	for(i = 0; i < size; i++)
		if (i != max_ind)
			sumDoubleKahan(sign[i] * max_sign * expm1(a[i]-max_double), &kexpSum); 
	sumDoubleKahan(max_sign * (2 * npos - size), &kexpSum);
	res = max_double+log(fabs(totalDoubleKahan(&kexpSum))); // absolute value here ???
	return res;
}
