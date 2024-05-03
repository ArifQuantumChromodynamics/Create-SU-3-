#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define SIZE 3

typedef struct{
	double real, imag;
}Complex;

Complex conjugate(Complex a){
	Complex result;
	result.real = a.real;
	result.imag = -a.imag;
	return result;
}

Complex multiplication(Complex a, Complex b){
	Complex result;
	result.real = a.real*b.real - a.imag*b.imag;
	result.imag = a.imag*b.real + a.real*b.imag;
	return result;
}

void SpecialUnitary3(Complex matrixGiven[SIZE][SIZE], Complex matrix_U[SIZE][SIZE]){
	Complex transpose[SIZE][SIZE];
	Complex X1[SIZE], X2[SIZE], X3[SIZE], Y1[SIZE], Y2[SIZE], Y3[SIZE];
	Complex prod_1 = {0, 0},  prod_2 = {0, 0},  prod_3 = {0, 0};
	Complex k_1 = {0, 0}, k_2 = {0, 0}, k_3 = {0, 0};
	float norm_1 = 0, norm_2 = 0, norm_3 = 0;
	Complex g = {0, 0}, h = {0, 0}, i = {0, 0};
	
	Complex orthoCheck1 = {0, 0}, orthoCheck2 = {0, 0}, orthoCheck3 = {0, 0}, normCheck1 = {0, 0}, normCheck2 = {0, 0}, normCheck3 = {0, 0}, detCheck = {0, 0};	

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			transpose[i][j] = matrixGiven[j][i];
		}
	}
	
	for (int j = 0; j <SIZE; j++){
		X1[j].real = transpose[0][j].real;
		X1[j].imag = transpose[0][j].imag;

		X2[j].real = transpose[1][j].real;
		X2[j].imag = transpose[1][j].imag;

		X3[j].real = transpose[2][j].real;
		X3[j].imag = transpose[2][j].imag;
	}

	printf("\nVectors that are the columns of the Given Matrix:\n");
	printf("X1 \t \tX2 \t  \tX3 \n");
	for (int k = 0; k < SIZE; k++){
		printf("%0.3f+(%0.3f)i\t%0.3f+(%0.3f)i\t%0.3f+(%0.3f)i\n", X1[k].real, X1[k].imag, X2[k].real, X2[k].imag, X3[k].real, X3[k].imag);
	}
	//Compute Y1
	for (int l = 0; l < SIZE; l++){
		Y1[l].real = X1[l].real;
		Y1[l].imag = X1[l].imag;
	}
	
	for (int m = 0; m < SIZE; m++){
		norm_1 += Y1[m].real*Y1[m].real + Y1[m].imag*Y1[m].imag;
		
		prod_1.real += multiplication(X2[m], conjugate(Y1[m])).real;
		prod_1.imag += multiplication(X2[m], conjugate(Y1[m])).imag;
	}
	
	k_1.real = prod_1.real/norm_1;
	k_1.imag = prod_1.imag/norm_1;

	//Compute Y2
	for(int m = 0; m < SIZE; m++){
		Y2[m].real = X2[m].real - multiplication(k_1, Y1[m]).real;
		Y2[m].imag = X2[m].imag - multiplication(k_1, Y1[m]).imag;	
	}
	
	for (int n = 0; n < SIZE; n++){
		norm_2 += Y2[n].real*Y2[n].real + Y2[n].imag*Y2[n].imag;
	
		prod_2.real += multiplication(X3[n], conjugate(Y1[n])).real;
		prod_2.imag += multiplication(X3[n], conjugate(Y1[n])).imag;

		prod_3.real += multiplication(X3[n], conjugate(Y2[n])).real;
		prod_3.imag += multiplication(X3[n], conjugate(Y2[n])).imag;
	}
	
	k_2.real = prod_2.real/norm_1;
	k_2.imag = prod_2.imag/norm_1;

	k_3.real = prod_3.real/norm_2;
	k_3.imag = prod_3.imag/norm_2;
	
	//Compute Y3
	for (int o = 0; o < SIZE; o++){
		Y3[o].real = X3[o].real - multiplication(k_2, Y1[o]).real - multiplication(k_3, Y2[o]).real;
		Y3[o].imag = X3[o].imag - multiplication(k_2, Y1[o]).imag - multiplication(k_3, Y2[o]).imag;
	}
	
	printf("\nThe Orthogonal Vectors obtained using Gram Schmidt process :- \n");
	printf("Y1: \t \t Y2: \t \t \t Y3: \n");
	for (int a = 0; a < SIZE; a++){
		printf("%0.3f+(%0.3f)i \t %0.3f+(%0.3f)i \t %0.3f+(%0.3f)i\n", Y1[a].real, Y1[a].imag, Y2[a].real, Y2[a].imag, Y3[a].real, Y3[a].imag);
	}
	
	for (int p = 0; p < SIZE; p++){
		norm_3 += Y3[p].real*Y3[p].real + Y3[p].imag*Y3[p].imag;
	}	

	for (int d = 0; d < SIZE; d++){
		Y1[d].real = Y1[d].real/sqrt(norm_1);
		Y1[d].imag = Y1[d].imag/sqrt(norm_1);
	
		Y2[d].real = Y2[d].real/sqrt(norm_2);
		Y2[d].imag = Y2[d].imag/sqrt(norm_2);

		Y3[d].real = Y3[d].real/sqrt(norm_3);
		Y3[d].imag = Y3[d].imag/sqrt(norm_3);
		
	}	
	printf("\nThe Orthogonal and Normalised Vectors are:- \n");
	printf("Y1 norm- \t Y2 norm- \t \t Y3 norm- \n");
	for (int b = 0; b < SIZE; b++){

		printf("%0.3f+(%0.3f)i \t %0.3f+(%0.3f)i \t %0.3f+(%0.3f)i\n", Y1[b].real, Y1[b].imag, Y2[b].real, Y2[b].imag, Y3[b].real, Y3[b].imag); 
	}
	
	//Set the determinant to 1
	
	g.real = multiplication(Y2[0], Y3[1]).real - multiplication(Y2[1], Y3[0]).real;
	g.imag = multiplication(Y2[0], Y3[1]).imag - multiplication(Y2[1], Y3[0]).imag;

	h.real = multiplication(Y1[1], Y3[0]).real - multiplication(Y1[0], Y3[1]).real;
	h.imag = multiplication(Y1[1], Y3[0]).imag - multiplication(Y1[0], Y3[1]).imag;

	i.real = multiplication(Y1[0], Y2[1]).real - multiplication(Y1[1], Y2[0]).real;
	i.imag = multiplication(Y1[0], Y2[1]).imag - multiplication(Y1[1], Y2[0]).imag;

	Y1[2].real = conjugate(g).real;
	Y1[2].imag = conjugate(g).imag;
	
	Y2[2].real = conjugate(h).real;
	Y2[2].imag = conjugate(h).imag;
	
	Y3[2].real = conjugate(i).real;
	Y3[2].imag = conjugate(i).imag;
	
	detCheck.real = multiplication(Y1[2], conjugate(Y1[2])).real + multiplication(Y2[2], conjugate(Y2[2])).real + multiplication(Y3[2], conjugate(Y3[2])).real;
	detCheck.imag = multiplication(Y1[2], conjugate(Y1[2])).imag + multiplication(Y2[2], conjugate(Y2[2])).imag + multiplication(Y3[2], conjugate(Y3[2])).imag;

	printf("\nCheck if the determinant is 1 of the resultant matrix:\n");
	printf("g*g + h*h + i*i = %0.3f+(%0.3f)i\n", detCheck.real, detCheck.imag);

	for(int t = 0; t < SIZE; t++){
		matrix_U[t][0] = Y1[t];
		matrix_U[t][1] = Y2[t];
		matrix_U[t][2] = Y3[t];
	}
	
	for (int b = 0; b < SIZE; b++){
		orthoCheck1.real += multiplication(Y1[b], conjugate(Y2[b])).real;
		orthoCheck1.imag += multiplication(Y1[b], conjugate(Y2[b])).imag;
		
		orthoCheck2.real += multiplication(Y2[b], conjugate(Y3[b])).real;
		orthoCheck2.imag += multiplication(Y2[b], conjugate(Y3[b])).imag;

		orthoCheck3.real += multiplication(Y3[b], conjugate(Y1[b])).real;
		orthoCheck3.imag += multiplication(Y3[b], conjugate(Y1[b])).imag;
	
		normCheck1.real += multiplication(Y1[b], conjugate(Y1[b])).real;
		normCheck1.imag += multiplication(Y1[b], conjugate(Y1[b])).imag;
		
		normCheck2.real += multiplication(Y2[b], conjugate(Y2[b])).real;
		normCheck2.imag += multiplication(Y2[b], conjugate(Y2[b])).imag;

		normCheck3.real += multiplication(Y3[b], conjugate(Y3[b])).real;
		normCheck3.imag += multiplication(Y3[b], conjugate(Y3[b])).imag;

	}
	
	printf("\nChecking the Orthogonality of the vectors-\n");
	printf("Y1 * bar(Y2): %0.03f+(%0.3f)i\n", orthoCheck1.real, orthoCheck1.imag);
	printf("Y2 * bar(Y3): %0.03f+(%0.3f)i\n", orthoCheck2.real, orthoCheck2.imag);
	printf("Y3 * bar(Y1): %0.03f+(%0.3f)i\n", orthoCheck3.real, orthoCheck3.imag);
	
	printf("\nChecking the Normalisation of the vectors-\n");
	printf("Y1 * bar(Y1): %0.03f+(%0.3f)i\n", normCheck1.real, normCheck1.imag);
	printf("Y2 * bar(Y2): %0.03f+(%0.3f)i\n", normCheck2.real, normCheck2.imag);
	printf("Y3 * bar(Y3): %0.03f+(%0.3f)i\n", normCheck3.real, normCheck3.imag);	
}
