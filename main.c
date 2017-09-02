#include "qustring.h"
#include <stdio.h>

unsigned char DeutschAlgorithm(const cmatnxm *u);

int main(int argc, char *argv[]){

	cmatnxm uf;
	cmatnxmInit(&uf, 4, 4);

	uf.m[0][0] = 1.0; uf.m[0][1] = 0.0; uf.m[0][2] = 0.0; uf.m[0][3] = 0.0;
	uf.m[1][0] = 0.0; uf.m[1][1] = 1.0; uf.m[1][2] = 0.0; uf.m[1][3] = 0.0;
	uf.m[2][0] = 0.0; uf.m[2][1] = 0.0; uf.m[2][2] = 1.0; uf.m[2][3] = 0.0;
	uf.m[3][0] = 0.0; uf.m[3][1] = 0.0; uf.m[3][2] = 0.0; uf.m[3][3] = 1.0;
	printf("Uf1:\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\nResult: %s\n\n", DeutschAlgorithm(&uf) == 0 ? "Constant" : "Balanced");

	uf.m[0][0] = 0.0; uf.m[0][1] = 1.0; uf.m[0][2] = 0.0; uf.m[0][3] = 0.0;
	uf.m[1][0] = 1.0; uf.m[1][1] = 0.0; uf.m[1][2] = 0.0; uf.m[1][3] = 0.0;
	uf.m[2][0] = 0.0; uf.m[2][1] = 0.0; uf.m[2][2] = 0.0; uf.m[2][3] = 1.0;
	uf.m[3][0] = 0.0; uf.m[3][1] = 0.0; uf.m[3][2] = 1.0; uf.m[3][3] = 0.0;
	printf("Uf2:\n0 1 0 0\n1 0 0 0\n0 0 0 1\n0 0 1 0\nResult: %s\n\n", DeutschAlgorithm(&uf) == 0 ? "Constant" : "Balanced");

	uf.m[0][0] = 1.0; uf.m[0][1] = 0.0; uf.m[0][2] = 0.0; uf.m[0][3] = 0.0;
	uf.m[1][0] = 0.0; uf.m[1][1] = 1.0; uf.m[1][2] = 0.0; uf.m[1][3] = 0.0;
	uf.m[2][0] = 0.0; uf.m[2][1] = 0.0; uf.m[2][2] = 0.0; uf.m[2][3] = 1.0;
	uf.m[3][0] = 0.0; uf.m[3][1] = 0.0; uf.m[3][2] = 1.0; uf.m[3][3] = 0.0;
	printf("Uf3:\n1 0 0 0\n0 1 0 0\n0 0 0 1\n0 0 1 0\nResult: %s\n\n", DeutschAlgorithm(&uf) == 0 ? "Constant" : "Balanced");

	uf.m[0][0] = 0.0; uf.m[0][1] = 1.0; uf.m[0][2] = 0.0; uf.m[0][3] = 0.0;
	uf.m[1][0] = 1.0; uf.m[1][1] = 0.0; uf.m[1][2] = 0.0; uf.m[1][3] = 0.0;
	uf.m[2][0] = 0.0; uf.m[2][1] = 0.0; uf.m[2][2] = 1.0; uf.m[2][3] = 0.0;
	uf.m[3][0] = 0.0; uf.m[3][1] = 0.0; uf.m[3][2] = 0.0; uf.m[3][3] = 1.0;
	printf("Uf4:\n0 1 0 0\n1 0 0 0\n0 0 1 0\n0 0 0 1\nResult: %s\n\n", DeutschAlgorithm(&uf) == 0 ? "Constant" : "Balanced");

	cmatnxmDelete(&uf);
	return 0;
}
