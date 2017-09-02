#include "qustring.h"

unsigned char DeutschAlgorithm(const cmatnxm *uf){

	// Create two qubits, |0> and |1>, and a temporary
	// qustring for holding certain function results.
	qustring q0, q1, t;
	qustringInitBasis(&q0, 1, 0);
	qustringInitBasis(&q1, 1, 1);

	// Create a new Hadamard gate.
	cmatnxm hadamard;
	cmatnxmNewHadamard(&hadamard);

	// Apply the Hadamard gate to both qubits.
	qustringMatrixMultKet(&hadamard, &q0, &t);
	qustringDelete(&q0); q0 = t;
	qustringMatrixMultKet(&hadamard, &q1, &t);
	qustringDelete(&q1); q1 = t;

	// Combine the two qubits into a single qustring.
	qustring q;
	qustringCombine(&q0, &q1, &q);
	qustringDelete(&q0);
	qustringDelete(&q1);

	// Apply the U gate to the qustring.
	qustringMatrixMultKet(uf, &q, &t);
	qustringDelete(&q);

	// Find the Kronecker product of the Hadamard gate
	// and an Identity matrix. The result is a gate which
	// we can multiply by the qustring in order to only
	// apply the Hadamard gate to the first qubit of the
	// qustring in the next operation.
	cmatnxm identity, hki;
	cmatnxmNewIdentity(&identity, 2);
	cmatnxmKroneckerProduct(&hadamard, &identity, &hki);
	cmatnxmDelete(&hadamard);
	cmatnxmDelete(&identity);

	// Apply the new gate to the qustring.
	qustringMatrixMultKet(&hki, &t, &q);
	cmatnxmDelete(&hki);
	qustringDelete(&t);

	// Measure the qustring.
	qustringMeasure(&q);

	// Get the value of the first qubit.
	unsigned char r = qustringGetQubitValue(&q, 0);

	// Clean up.
	qustringDelete(&q);

	// Return.
	return r;

}
