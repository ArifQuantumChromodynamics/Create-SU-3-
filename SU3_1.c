#include "SU3_0.c"

void checkSpecialUnitary3(Complex matrix_U[SIZE][SIZE]){
	Complex conjTrap[SIZE][SIZE], result[SIZE][SIZE];

	Complex det = {0, 0}, aei = {0, 0}, bfg = {0, 0}, cdh = {0, 0};
	Complex afh = {0, 0}, bdi = {0, 0}, ceg = {0, 0};
 
	//Create the Conjugate Transpose Matrix

	for (int e = 0; e < SIZE; e++){
		for (int f = 0; f < SIZE; f++){
			conjTrap[e][f] = conjugate(matrix_U[f][e]);
		}
	}
	printf("\nThe conjugate transpose Matrix is—\n");
	for (int a = 0; a < SIZE; a++){
		for (int b = 0; b < SIZE; b++){
			printf("%0.3f + (%0.3f)i\t", conjTrap[a][b].real, conjTrap[a][b].imag);
		}
	printf("\n");
	}

	//Calculate the product of U^{dagg} and U

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			result[i][j].real = 0;
			result[i][j].imag = 0;
			for (int k = 0; k < SIZE; k++) {
				Complex term = multiplication(matrix_U[k][j], conjTrap[i][k]);
				result[i][j].real += term.real;
				result[i][j].imag += term.imag;
			}
		}
	}
	
	printf("\nProduct of the matrix SU(3) and it's conjugate transpose:\n");
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			printf("%0.3f+(%0.3f)i   ", result[i][j].real, result[i][j].imag);
        }
        printf("\n");
	}

	printf("\nDeterminant of the SU(3) matrix: ");
	aei = multiplication(matrix_U[0][0], multiplication(matrix_U[1][1], matrix_U[2][2]));
	bfg = multiplication(matrix_U[0][1], multiplication(matrix_U[1][2], matrix_U[2][0]));
	cdh = multiplication(matrix_U[0][2], multiplication(matrix_U[1][0], matrix_U[2][1]));

	afh = multiplication(matrix_U[0][0], multiplication(matrix_U[1][2], matrix_U[2][1]));
    	bdi = multiplication(matrix_U[0][1], multiplication(matrix_U[1][0], matrix_U[2][2]));
	ceg = multiplication(matrix_U[0][2], multiplication(matrix_U[1][1], matrix_U[2][0]));

	 det.real = aei.real + bfg.real + cdh.real - (afh.real + bdi.real + ceg.real);
    det.imag = aei.imag + bfg.imag + cdh.imag - (afh.imag + bdi.imag + ceg.imag);
	printf("%0.3f + (%0.3f)i\n",  det.real,  det.imag);
}