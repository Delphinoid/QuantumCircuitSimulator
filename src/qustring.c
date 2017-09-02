#include "qustring.h"
#include <math.h>
#include <float.h>
#include <string.h>

double complex qubitStateProbability(const double complex a){
	/*
	** The probability of a state is |a|^2.
	*/
	double complex sqcabs = cabs(a);
	return sqcabs*sqcabs;
}

unsigned char qubitIsValid(const qubit *b){
	return fabs(qubitStateProbability((*b)[0]) + qubitStateProbability((*b)[1]))
	       > 1.0 - DBL_EPSILON;
}

void qubitBlochSphere(const qubit *b, float pitch, float yaw, float roll){
	/*
	** Finds the corresponding bloch sphere representation of
	** the specified qubit.
	*/
}

unsigned char qustringInit(qustring *q, const size_t n){
	/*
	** Initializes an n-qubit qustring.
	*/
	q->n = n<<1;
	q->v = malloc(q->n*sizeof(double complex));
	if(q->v == NULL){
		q->n = 0;
		return 0;
	}
	return 1;
}

unsigned char qustringInitBasis(qustring *q, const size_t n, const size_t s){
	/*
	** Initializes an n-qubit qustring in its sth basis state.
	*/
	q->n = n<<1;
	q->v = malloc(q->n*sizeof(double complex));
	if(q->v == NULL){
		q->n = 0;
		return 0;
	}
	if(s < q->n){
		qustringCollapse(q, s);
	}
	return 1;
}

unsigned char qustringAddControl(qustring *q, const qubit *b){
	/*
	** Appends a control qubit to a qustring. Used for certain
	** gates that require multiple parameters.
	*/
	double complex *tempBuffer = realloc((void *)q->v, (q->n+2)*sizeof(double complex));
	if(tempBuffer == NULL){
		return 0;
	}
	q->v = tempBuffer;
	q->v[q->n]   = (*b)[0];
	q->v[++q->n] = (*b)[1];
	++q->n;
	return 1;
}

unsigned char qustringIsValid(const qustring *q){
	/*
	** Make sure the unitarity property for the given
	** qustring holds (the probabilities for each possible
	** state add up to 1). That is, evaluate the following
	** condition:
	**
	** q->n
	**  __
	**  \ '
	**  /_, |q[i]|^2 = 1
	**
	**  i=0
	**
	** (Yes, I DO like drawing equations in ASCII!)
	*/
	double sum = 0.0;
	size_t i;
	for(i = 0; i < q->n; ++i){
		sum += qubitStateProbability(q->v[i]);
	}
	return fabs(sum) > 1.0 - DBL_EPSILON;
}

unsigned char qustringGetQubit(const qustring *q, size_t n, qubit *r){
	n = n<<1;
	if(++n < q->n){
		(*r)[1] = q->v[n];
		(*r)[0] = q->v[--n];
		return 1;
	}
	return 0;
}

unsigned char qustringGetQubitValue(const qustring *q, size_t n){
	/*
	** The result of this function is only
	** meaningful on a measured qubit.
	*/
	n = n<<1;
	if(++n < q->n){
		return (unsigned char)creal(qubitStateProbability(q->v[n]));
	}
	return 0;
}

void qustringCollapse(qustring *q, const size_t s){
	/*
	** Collapses a qustring into a pure (basis) state
	** representing s.
	*/
	memset((void *)q->v, 0, q->n*sizeof(q->v[0]));
	q->v[s] = 1.0;
}

unsigned int qustringMeasure(qustring *q){
	/*
	** Measure the state of a qustring. For each individual
	** state, |state|^2 is its probability. By converting
	** the probability of each state to a percentage of
	** RAND_MAX and removing it from the result of rand(),
	** we can approximate such an effect.
	*/
	// Using a random float / double instead of an int makes
	// this a lot easier.
	double d = rand() / (RAND_MAX + 1.0);
	size_t i;
	for(i = 0; i < q->n; ++i){
		d -= qubitStateProbability(q->v[i]) * RAND_MAX;
		if(d > 0.0){
			q->v[i] = 0.0;
		}else{
			// We have evaluated a state for the qubit,
			// collapse it into a pure state and return
			// the result.
			memset((void *)(q->v+i+1), 0, (q->n-i-1)*sizeof(q->v[0]));
			q->v[i] = 1.0;
			return i;
		}
	}
	// Function failed, qustring is most likely invalid or
	// there was a precision error. Assume the latter.
	return 0;
}

unsigned char qustringCombine(const qustring *q1, const qustring *q2, qustring *r){
	/*
	** The combination of two qustrings is equivalent to
	** their Tensor product.
	*/
	if(!qustringInit(r, (q1->n*q2->n)>>1)){
		return 0;
	}
	size_t i, j;
	for(i = 0; i < q1->n; ++i){
		for(j = 0; j < q2->n; ++j){
			r->v[(i * q2->n) + j] = q1->v[i] * q2->v[j];
		}
	}
	return 1;
}

unsigned char qustringMatrixMultBra(const qustring *q, const cmatnxm *m, qustring *r){
	/*
	** Multiplies a bra-vector by a matrix (such as a
	** quantum logic gate).
	*/
	if(q->n == m->r){
		if(!qustringInit(r, m->c>>1)){
			return 0;
		}
		size_t i, j;
		for(i = 0; i < r->n; ++i){
			r->v[i] = 0.0;
			for(j = 0; j < q->n; ++j){
				r->v[i] += q->v[j] * m->m[j][i];
			}
		}
	}
	return 1;
}

unsigned char qustringMatrixMultKet(const cmatnxm *m, const qustring *q, qustring *r){
	/*
	** Multiplies a ket-vector by a matrix (such as a
	** quantum logic gate).
	*/
	if(m->c == q->n){
		if(!qustringInit(r, m->r>>1)){
			return 0;
		}
		size_t i, j;
		for(i = 0; i < r->n; ++i){
			r->v[i] = 0.0;
			for(j = 0; j < q->n; ++j){
				r->v[i] += m->m[i][j] * q->v[j];
			}
		}
	}
	return 1;
}

void qustringDelete(qustring *q){
	if(q->v != NULL){
		free(q->v);
	}
}
