#include "cmatnxm.h"
#include <math.h>
#include <string.h>

#define SQRT_HALF   0.7071067811865476
#define SQRT_HALF_I 0.7071067811865476*I

unsigned char cmatnxmInit(cmatnxm *m, const size_t r, const size_t c){
	double complex **mat = malloc(r * sizeof(double complex *));
	if(mat == NULL){
		return 0;
	}
	size_t i;
	for(i = 0; i < r; ++i){
		mat[i] = malloc(c * sizeof(double complex));
		if(mat[i] == NULL){
			/*
			** If malloc() failed, free each previous
			** element and return.
			*/
			size_t j;
			for(j = i; j > 0; --i){
				free(mat[j-1]);
			}
			free(mat);
			return 0;
		}
	}
	m->r = r;
	m->c = c;
	m->m = mat;
	return 1;
}

void cmatnxmIdentity(cmatnxm *m){
	size_t i, j;
	for(i = 0; i < m->r; ++i){
		for(j = 0; j < m->c; ++j){
			if(i == j){
				m->m[i][j] = 1.0;
			}else{
				m->m[i][j] = 0.0;
			}
		}
	}
}

unsigned char cmatnxmNewIdentity(cmatnxm *m, const size_t n){
	if(!cmatnxmInit(m, n, n)){
		return 0;
	}
	cmatnxmIdentity(m);
	return 1;
}

unsigned char cmatnxmNewHadamard(cmatnxm *m){
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = SQRT_HALF; m->m[0][1] =  SQRT_HALF;
	m->m[1][0] = SQRT_HALF; m->m[1][1] = -SQRT_HALF;
	return 1;
}

unsigned char cmatnxmNewPhaseShift(cmatnxm *m, const double complex s){
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = cexp(I * s);
	return 1;
}

unsigned char cmatnxmNewPauliSpinX(cmatnxm *m){
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 0.0; m->m[0][1] = 1.0;
	m->m[1][0] = 1.0; m->m[1][1] = 0.0;
	return 1;
}

unsigned char cmatnxmNewPauliSpinY(cmatnxm *m){
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 0.0; m->m[0][1] =   -I;
	m->m[1][0] =   I; m->m[1][1] =  0.0;
	return 1;
}

unsigned char cmatnxmNewPauliSpinZ(cmatnxm *m){
	// return cmatnxmNewPhaseShift(m, M_PI);
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] =  0.0;
	m->m[1][0] = 0.0; m->m[1][1] = -1.0;
	return 1;
}

unsigned char cmatnxmNewPhaseI(cmatnxm *m){
	// return cmatnxmNewPhaseShift(m, M_PI_2);
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] =   I;
	return 1;
}

unsigned char cmatnxmNewPiOver8(cmatnxm *m){
	//return cmatnxmNewPhaseShift(m, M_PI_4);
	if(!cmatnxmInit(m, 2, 2)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = SQRT_HALF + SQRT_HALF_I;
	return 1;
}

unsigned char cmatnxmNewControlledNot(cmatnxm *m){
	if(!cmatnxmInit(m, 4, 4)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0; m->m[0][3] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0; m->m[1][3] = 0.0;
	m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 0.0; m->m[2][3] = 1.0;
	m->m[3][0] = 0.0; m->m[3][1] = 0.0; m->m[3][2] = 1.0; m->m[3][3] = 0.0;
	return 1;
}

unsigned char cmatnxmNewControlledZ(cmatnxm *m){
	if(!cmatnxmInit(m, 4, 4)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0; m->m[0][3] =  0.0;
	m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0; m->m[1][3] =  0.0;
	m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 1.0; m->m[2][3] =  0.0;
	m->m[3][0] = 0.0; m->m[3][1] = 0.0; m->m[3][2] = 0.0; m->m[3][3] = -1.0;
	return 1;
}

unsigned char cmatnxmNewControlledShift(cmatnxm *m){
	if(!cmatnxmInit(m, 4, 4)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0; m->m[0][3] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0; m->m[1][3] = 0.0;
	m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 1.0; m->m[2][3] = 0.0;
	m->m[3][0] = 0.0; m->m[3][1] = 0.0; m->m[3][2] = 0.0; m->m[3][3] = I;
	return 1;
}

unsigned char cmatnxmNewToffoli(cmatnxm *m){
	if(!cmatnxmInit(m, 8, 8)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0; m->m[0][3] = 0.0; m->m[0][4] = 0.0; m->m[0][5] = 0.0; m->m[0][6] = 0.0; m->m[0][7] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0; m->m[1][3] = 0.0; m->m[1][4] = 0.0; m->m[1][5] = 0.0; m->m[1][6] = 0.0; m->m[1][7] = 0.0;
	m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 1.0; m->m[2][3] = 0.0; m->m[2][4] = 0.0; m->m[2][5] = 0.0; m->m[2][6] = 0.0; m->m[2][7] = 0.0;
	m->m[3][0] = 0.0; m->m[3][1] = 0.0; m->m[3][2] = 0.0; m->m[3][3] = 1.0; m->m[3][4] = 0.0; m->m[3][5] = 0.0; m->m[3][6] = 0.0; m->m[3][7] = 0.0;
	m->m[4][0] = 0.0; m->m[4][1] = 0.0; m->m[4][2] = 0.0; m->m[4][3] = 0.0; m->m[4][4] = 1.0; m->m[4][5] = 0.0; m->m[4][6] = 0.0; m->m[4][7] = 0.0;
	m->m[5][0] = 0.0; m->m[5][1] = 0.0; m->m[5][2] = 0.0; m->m[5][3] = 0.0; m->m[5][4] = 0.0; m->m[5][5] = 1.0; m->m[5][6] = 0.0; m->m[5][7] = 0.0;
	m->m[6][0] = 0.0; m->m[6][1] = 0.0; m->m[6][2] = 0.0; m->m[6][3] = 0.0; m->m[6][4] = 0.0; m->m[6][5] = 0.0; m->m[6][6] = 0.0; m->m[6][7] = 1.0;
	m->m[7][0] = 0.0; m->m[7][1] = 0.0; m->m[7][2] = 0.0; m->m[7][3] = 0.0; m->m[7][4] = 0.0; m->m[7][5] = 0.0; m->m[7][6] = 1.0; m->m[7][7] = 0.0;
	return 1;
}

unsigned char cmatnxmNewFredkin(cmatnxm *m){
	if(!cmatnxmInit(m, 8, 8)){
		return 0;
	}
	m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0; m->m[0][3] = 0.0; m->m[0][4] = 0.0; m->m[0][5] = 0.0; m->m[0][6] = 0.0; m->m[0][7] = 0.0;
	m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0; m->m[1][3] = 0.0; m->m[1][4] = 0.0; m->m[1][5] = 0.0; m->m[1][6] = 0.0; m->m[1][7] = 0.0;
	m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 1.0; m->m[2][3] = 0.0; m->m[2][4] = 0.0; m->m[2][5] = 0.0; m->m[2][6] = 0.0; m->m[2][7] = 0.0;
	m->m[3][0] = 0.0; m->m[3][1] = 0.0; m->m[3][2] = 0.0; m->m[3][3] = 1.0; m->m[3][4] = 0.0; m->m[3][5] = 0.0; m->m[3][6] = 0.0; m->m[3][7] = 0.0;
	m->m[4][0] = 0.0; m->m[4][1] = 0.0; m->m[4][2] = 0.0; m->m[4][3] = 0.0; m->m[4][4] = 1.0; m->m[4][5] = 0.0; m->m[4][6] = 0.0; m->m[4][7] = 0.0;
	m->m[5][0] = 0.0; m->m[5][1] = 0.0; m->m[5][2] = 0.0; m->m[5][3] = 0.0; m->m[5][4] = 0.0; m->m[5][5] = 0.0; m->m[5][6] = 1.0; m->m[5][7] = 0.0;
	m->m[6][0] = 0.0; m->m[6][1] = 0.0; m->m[6][2] = 0.0; m->m[6][3] = 0.0; m->m[6][4] = 0.0; m->m[6][5] = 1.0; m->m[6][6] = 0.0; m->m[6][7] = 0.0;
	m->m[7][0] = 0.0; m->m[7][1] = 0.0; m->m[7][2] = 0.0; m->m[7][3] = 0.0; m->m[7][4] = 0.0; m->m[7][5] = 0.0; m->m[7][6] = 0.0; m->m[7][7] = 1.0;
	return 1;
}

unsigned char cmatnxmCopy(const cmatnxm *m1, cmatnxm *m2){
	if(!cmatnxmInit(m2, m1->r, m1->c)){
		return 0;
	}
	size_t i;
	for(i = 0; i < m1->r; ++i){
		m2->m[i][0] = m1->m[i][0];
		memcpy(m2->m[i], m1->m[i], m1->c*sizeof(m2->m[i][0]));
	}
	return 1;
}

unsigned char cmatnxmMMultM(const cmatnxm *m1, const cmatnxm *m2, cmatnxm *r){
	if(m1->c == m2->r){
		if(!cmatnxmInit(r, m1->r, m2->c)){
			return 0;
		}
		size_t i, j, k;
		for(i = 0; i < m1->r; ++i){
			for(j = 0; j < m2->c; ++j){
				r->m[i][j] = 0.0;
				for(k = 0; k < m1->c; ++k){
					r->m[i][j] += m1->m[i][k] * m2->m[k][j];
				}
			}
		}
	}
	return 1;
}

unsigned char cmatnxmKroneckerProduct(const cmatnxm *m1, const cmatnxm *m2, cmatnxm *r){
	if(!cmatnxmInit(r, m1->r*m2->r, m1->c*m2->c)){
		return 0;
	}
	size_t i, j, k, l;
	for(i = 0; i < m1->r; ++i){
		for(j = 0; j < m1->c; ++j){
			for(k = 0; k < m2->r; ++k){
				for(l = 0; l < m2->c; ++l){
					r->m[i*m2->r+k][j*m2->c+j] = m1->m[i][j] * m2->m[k][l];
				}
			}
		}
	}
	return 1;
}

unsigned char cmatnxmKroneckerProductExponential(const cmatnxm *m, unsigned int n, cmatnxm *r){
	*r = *m;
	for(; n > 0; --n){
		cmatnxm temp;
		if(!cmatnxmCopy(r, &temp)){
			cmatnxmDelete(r);
			return 0;
		}
		cmatnxmDelete(r);
		if(!cmatnxmKroneckerProduct(m, &temp, r)){
			cmatnxmDelete(&temp);
			return 0;
		}
		cmatnxmDelete(&temp);
	}
	return 1;
}

void cmatnxmDelete(cmatnxm *m){
	if(m->m != NULL){
		size_t i;
		for(i = 0; i < m->r; ++i){
			free(m->m[i]);
		}
		free(m->m);
	}
	m->r = 0;
	m->c = 0;
}
