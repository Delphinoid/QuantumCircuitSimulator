#ifndef QUSTRING_H
#define QUSTRING_H

#include "cmatnxm.h"

/*
** A qubit can store the basis states |0> or |1>, or a
** superposition of both. These states are represented
** by complex numbers, where the following holds:
**
**  n
** __
** \ '
** /_, |q[i]|^2 = 1
**
** i=0
**
** These numbers represent the probability amplitudes of
** either state. That is, |q[0]|^2 is the probability that
** the qubit is in an unset (0) state, while |q[1]|^2 is
** the probability that it is in a set (1) state.
*/
typedef double complex qubit[2];

double complex qubitStateProbability(const double complex a);
unsigned char qubitIsValid(const qubit *b);
void qubitBlochSphere(const qubit *b, float pitch, float yaw, float roll);


/*
** The qustring (which can, more specifically, be described
** as a qubit register) represents a qubit string, similar
** to how a classical binary string represents a string of
** classical binary bits. We store them as a one-dimensional
** vector, which may be interpreted as a bra-vector (row
** vector) or a ket-vector (column vector) depending on the
** functions used. The probability amplitudes of the states
** for each bit are stored contiguously in order to make
** gates evaluation easier.
*/
typedef struct {
	size_t n;           // Vector size (number of qubits * 2, or the number of possible states)
	double complex *v;  // Qubit state vector
} qustring;

unsigned char qustringInit(qustring *q, const size_t n);
unsigned char qustringInitBasis(qustring *q, const size_t n, const size_t s);
unsigned char qustringAddControl(qustring *q, const qubit *b);
unsigned char qustringIsValid(const qustring *q);
unsigned char qustringGetQubit(const qustring *q, size_t n, qubit *r);
unsigned char qustringGetQubitValue(const qustring *q, size_t n);
void qustringCollapse(qustring *q, const size_t s);
unsigned int qustringMeasure(qustring *q);
unsigned char qustringCombine(const qustring *q1, const qustring *q2, qustring *r);
unsigned char qustringMatrixMultBra(const qustring *q, const cmatnxm *m, qustring *r);
unsigned char qustringMatrixMultKet(const cmatnxm *m, const qustring *q, qustring *r);
void qustringDelete(qustring *q);

#endif
