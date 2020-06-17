#include <iostream>
#include <iomanip>
#include <fstream>


using namespace std;
const int N = 5, M = 5;
ifstream input;
ofstream output;


class Matrix {
public:
	void readmat(double A[M][N])
	{
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				input >> A[i][j];
			}
		}
	}

	void copymat(double A[M][N], double D[M][N])
	{
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				D[i][j] = A[i][j];
			}
		}
	}

	void identmat(double B[M][N])
	{

		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i == j)
				{
					B[i][j] = 1;
				}
				else
				{
					B[i][j] = 0;
				}
			}
		}
	}

	void multmat(double D[M][N], double A[M][N], double C[M][N])
	{


		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				double sum = 0;
				for (int k = 0; k < M; k++)
				{
					C[i][j] = A[i][k] * D[k][j];
					sum = sum + C[i][j];

				}
				output << sum << " ";
			}
			output << endl;
		}
	}

	void printmat(double A[M][N])
	{
		double D[M][N];

		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				output << A[i][j] << " ";
			}
			output << "\n";
		}
	}


	void invMat(double A[M][N], double B[N][N]) {
		if (M == N) {
			for (int i = 0; i < N; i++) {
				for (int k = i; k < N; k++) {
					if (A[k][i] == 0) {					//if A[k][i] is equal to zero, swap so that it can be used to decrement (cannot use zero as a decrementing factor)
						for (int j = 0; j < N; j++) {		//do the same with matrix B, since any matrix operation on A must also be done on B
							float temp = A[i][j];
							float temp2 = B[i][j];
							A[i][j] = A[k][j];
							B[i][j] = B[k][j];

							A[k][j] = temp;
							B[k][j] = temp2;
						}
					}
				}
				double val = A[i][i];
				for (int c = 0; c < N; c++) {
					A[i][c] = A[i][c] / val;			//divide by the diagonal element in each row ( (1,1), (2,2), etc.)
					B[i][c] = B[i][c] / val;			//same operation on B
				}
				for (int z = 0; z < N; z++) {
					if (i != z) {									// for all values that are not diagonal matrix elements, subtract 
						double val2 = A[z][i];						// A[i][c] (remember that the whole row was divided by the corresponding diagonal element)
						for (int c = 0; c < N; c++) {			    // times  val2 so that the resulting value is zero in every non-diagonal element, yielding the
							A[z][c] = A[z][c] - (A[i][c] * val2);	// identity matrix in A. Apply the same formulas to B.
							B[z][c] = B[z][c] - (B[i][c] * val2);
						}
					}
				}

			}
		}
	}


};
int main(){
input.open("input1.txt");
output.open("output1.txt");
output << fixed << showpoint << setprecision(3);

Matrix myMatrix;

double A[M][N];
myMatrix.readmat(A);

output << "A:" << endl;
myMatrix.printmat(A);
output << endl;

double D[M][N];
myMatrix.copymat(A, D);

double B[N][N];
myMatrix.identmat(B);

output << "B:" << endl;
myMatrix.invMat(A, B);
myMatrix.printmat(B);
output << endl;



double C[N][N];
output << "C (multmat): " << endl;
myMatrix.multmat(D, B, C);
output << endl;


}






