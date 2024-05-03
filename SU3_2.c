#include "SU3_1.c"

int main(){
	Complex matrixGiven[SIZE][SIZE];
	
	srand(time(NULL));
	for (int i = 0; i < SIZE; i++){
		for (int j = 0; j < SIZE; j++){
			matrixGiven[i][j].real = (double)rand()/RAND_MAX;
			matrixGiven[i][j].imag = (double)rand()/RAND_MAX;
		}
	}
	printf("\nThe Given Matrix: \n");
	for (int a = 0; a < SIZE; a++){
		for (int b = 0; b < SIZE; b++){
			printf("%0.3f + (%0.3f)i \t", matrixGiven[a][b].real, matrixGiven[a][b].imag);
		}
	printf("\n");
	}
	
	Complex matrix_U[SIZE][SIZE];
	SpecialUnitary3(matrixGiven, matrix_U);
	
	printf("\nThe SU(3) Matrix is—\n");
	for (int a = 0; a < SIZE; a++){
		for (int b = 0; b < SIZE; b++){
			printf("%0.3f + (%0.3f)i\t", matrix_U[a][b].real, matrix_U[a][b].imag);
		}
	printf("\n");
	}
	
	checkSpecialUnitary3(matrix_U);
return 0;
}