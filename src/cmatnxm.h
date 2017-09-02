#ifndef CMATNXM_H
#define CMATNXM_H

#include <complex.h>
#include <stdlib.h>

typedef struct {
	size_t r, c;
	double complex **m;
} cmatnxm;

unsigned char cmatnxmInit(cmatnxm *m, const size_t r, const size_t c);
void cmatnxmIdentity(cmatnxm *m);
unsigned char cmatnxmNewIdentity(cmatnxm *m, const size_t n);
unsigned char cmatnxmNewHadamard(cmatnxm *m);
unsigned char cmatnxmNewPhaseShift(cmatnxm *m, const double complex s);
unsigned char cmatnxmNewPauliSpinX(cmatnxm *m);
unsigned char cmatnxmNewPauliSpinY(cmatnxm *m);
unsigned char cmatnxmNewPauliSpinZ(cmatnxm *m);
unsigned char cmatnxmNewPhaseI(cmatnxm *m);
unsigned char cmatnxmNewPiOver8(cmatnxm *m);
unsigned char cmatnxmNewControlledNot(cmatnxm *m);
unsigned char cmatnxmNewControlledZ(cmatnxm *m);
unsigned char cmatnxmNewControlledShift(cmatnxm *m);
unsigned char cmatnxmNewToffoli(cmatnxm *m);
unsigned char cmatnxmNewFredkin(cmatnxm *m);
unsigned char cmatnxmKroneckerProduct(const cmatnxm *m1, const cmatnxm *m2, cmatnxm *r);
unsigned char cmatnxmKroneckerProductExponential(const cmatnxm *m, unsigned int n, cmatnxm *r);
unsigned char cmatnxmMMultM(const cmatnxm *m1, const cmatnxm *m2, cmatnxm *r);
void cmatnxmDelete(cmatnxm *m);

#endif
