#ifndef KahanF
#define KahanF


#define sign(x) ((x)/fabs((x)))

typedef struct COMP_DOUBLE_KAHAN {
	double sum, c, cc, cs, ccs;
} TypeDoubleKahan;

#ifdef __cplusplus
extern "C" {
#endif

void initDoubleKahan(TypeDoubleKahan *ka);
void sumDoubleKahan(double x, TypeDoubleKahan *ka);
double totalDoubleKahan(TypeDoubleKahan *ka);
double sumSignedLogTableKahan(double a[], int sign[], int size);

#ifdef __cplusplus
}
#endif

#endif
