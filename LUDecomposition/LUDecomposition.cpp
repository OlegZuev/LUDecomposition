#include <iostream>
#include <fstream>
#include <cmath>
#include "LUDecomposition.h"
#include <string>
#include <iomanip>

using namespace std;

const int PRECISION = 7;
const double EPS = 1E-6;

int main() {
	double** matrix_A = nullptr;
	double** matrix_PA = nullptr;
	double** matrix_inverseA = nullptr;
	double** matrix_U = nullptr;
	double** matrix_L = nullptr;
	string variant = "8b";

	int n;
	ifstream fin("../" + variant + ".txt");
	ofstream fout("../" + variant + "_output" + ".txt");
	fin >> n;

	ostream& out = fout;

	int* indexes = new int[n];
	for (int i = 0; i < n; ++i) {
		indexes[i] = i;
	}

	double true_x[] = { 1, 2, 3, 4 };
	double* b = new double[4]{ 0 };

	allocate_matrix(matrix_A, n);
	allocate_matrix(matrix_PA, n);
	allocate_matrix(matrix_inverseA, n);
	allocate_matrix(matrix_U, n);
	allocate_matrix(matrix_L, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> matrix_A[i][j];
		}
	}

	fin.close();

	copy_matrix(matrix_U, matrix_A, n);
	copy_matrix(matrix_PA, matrix_A, n);

	out << fixed << setprecision(7) << "Variant=" << variant << endl;
	out << "x:" << endl;
	print_array(true_x, n, out);
	out << "A:" << endl;
	print_matrix(matrix_PA, n, out);
	out << "b:" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			b[i] += matrix_A[i][j] * true_x[j];
		}
	}
	print_array(b, n, out);

	int permutation_count;
	find_LU(matrix_PA, matrix_U, matrix_L, indexes, n, out, permutation_count);

	for (int i = 0; i < n; i++) {
		out << indexes[i] + 1 << " ";
	}
	out << endl;

	double** tmp = multiply_matrix(matrix_L, matrix_U, n);
	double** error = subtraction_matrix(tmp, matrix_PA, n);
	delete_matrix(tmp, n);

	out << "L*U - P*A" << endl;
	print_matrix(error, n, out << scientific << setprecision(PRECISION * 2));
	out << fixed << setprecision(PRECISION);

	double det_A = 1;
	for (int i = 0; i < n; i++) {
		det_A *= matrix_L[i][i];
	}

	det_A *= permutation_count % 2 == 0 ? 1 : -1;
	out << "Determinant=" << det_A << endl;
	out << "x:" << endl;
	print_array(true_x, n, out);
	out << "A:" << endl;
	print_matrix(matrix_A, n, out);

	find_inverse_matrix(matrix_inverseA, matrix_U, matrix_L, n, indexes);

	out << "A^(-1):" << endl;
	print_matrix(matrix_inverseA, n, out);

	out << "A*A^(-1)" << endl;
	tmp = multiply_matrix(matrix_A, matrix_inverseA, n);
	print_matrix(tmp, n, out << scientific << setprecision(PRECISION * 2));
	out << fixed << setprecision(PRECISION);

	out << "The condition number of the matrix A=" << endl;
	out << "Cond-1: " << find_n1(matrix_A, n) * find_n1(matrix_inverseA, n) << endl;
	out << "Cond-2: " << find_n2(matrix_A, n) * find_n2(matrix_inverseA, n) << endl;
	out << "Cond-3: " << find_n3(matrix_A, n) * find_n3(matrix_inverseA, n) << endl;

	double result_x[4]{ 0 };

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result_x[i] += matrix_inverseA[i][j] * b[j];
		}
	}

	out << "A*x-b:" << endl;
	double some[4]{ 0 };
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			some[i] += matrix_A[i][j] * result_x[j];
		}
		out << scientific << setprecision(PRECISION * 2) << some[i] - b[i] << " ";
	}

	out << endl << "Error: " << endl;
	for (int i = 0; i < n; ++i) {		
		out << scientific << setprecision(PRECISION * 2) << true_x[i] - result_x[i] << " ";
	}
	out << endl;

	delete_matrix(tmp, n);
	delete_matrix(matrix_PA, n);
	delete_matrix(matrix_inverseA, n);
	delete_matrix(matrix_A, n);
	delete_matrix(matrix_U, n);
	delete_matrix(matrix_L, n);
	delete_matrix(error, n);
	delete[] indexes;
	delete[] b;
	fout.close();
}

void allocate_matrix(double**& matrix, int& n) {
	matrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new double[n] {0};
	}
}

void delete_matrix(double** matrix, int& n) {
	for (int i = 0; i < n; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

void copy_matrix(double** matrix_to, double** matrix_from, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			matrix_to[i][j] = matrix_from[i][j];
		}
	}
}

void print_matrix(double** matrix, int n, ostream& ostr) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			ostr << matrix[i][j] << " ";
		}
		ostr << endl;
	}
	ostr << endl;
}

void print_array(double* arr, int n, ostream& ostr) {
	for (int i = 0; i < n; ++i) {
		ostr << arr[i] << " ";
	}
	ostr << endl;
}


double** multiply_matrix(double** matrix_a, double** matrix_b, int n) {
	double** new_matrix = nullptr;
	allocate_matrix(new_matrix, n);
	for (int i = 0; i < n; ++i) {
		for (int k = 0; k < n; ++k) {
			for (int j = 0; j < n; ++j) {
				new_matrix[i][k] += matrix_a[i][j] * matrix_b[j][k];
			}
		}
	}

	return new_matrix;
}

double** subtraction_matrix(double** matrix_a, double** matrix_b, int n) {
	double** result = nullptr;
	allocate_matrix(result, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = matrix_a[i][j] - matrix_b[i][j];
		}
	}

	return result;
}

void find_LU(double** matrix_PA, double** matrix_U, double** matrix_L, int* indexes, int n, ostream& out, int& permutation_count) {
	permutation_count = 0;
	int rank = n;
	for (int i = 0; i < n; ++i) {
		double column_max = matrix_U[i][i];
		int k = i;
		for (int j = 0; j < n; ++j) {
			if (abs(column_max) < abs(matrix_U[j][i])) {
				column_max = matrix_U[j][i];
				k = j;
			}
		}

		if (fabs(matrix_PA[k][i]) < EPS) {
			rank--;
		}

		if (i != k) {
			permutation_count++;
		}

		swap(matrix_U[i], matrix_U[k]);
		swap(indexes[i], indexes[k]);
		swap(matrix_PA[i], matrix_PA[k]);

		double multiplier = matrix_U[i][i];
		for (int j = 0; j < n; ++j) {
			matrix_U[i][j] = matrix_U[i][j] / multiplier;
		}

		for (int b = i + 1; b < n; ++b) {
			multiplier = matrix_U[b][i];
			for (int j = i; j < n; ++j) {
				matrix_U[b][j] = matrix_U[b][j] - multiplier * matrix_U[i][j];
			}
		}

		//out << "permutation count= " << permutation_count << " m= " << k + 1 << endl;
		out << "k= " << i + 1 << " m= " << k + 1 << endl;
		out << "U:" << endl;
		print_matrix(matrix_U, n, out);

		for (int j = 0; j <= i; ++j) {
			double vector_mult_result = 0;
			for (int k = 0; k < i; ++k) {
				vector_mult_result += matrix_L[i][k] * matrix_U[k][j];
			}

			matrix_L[i][j] = matrix_PA[i][j] - vector_mult_result;
		}

		out << "L:" << endl;
		print_matrix(matrix_L, n, out);
	}

	out << "Rank= " << rank << endl;
}

void find_inverse_matrix(double** matrix_inverseA, double** matrix_U, double** matrix_L, int n, int* indexes) {
	double* y = new double[n] {0};
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double sum = 0;
			for (int k = 0; k < j; ++k) {
				sum += matrix_L[j][k] * y[k];
			}

			y[j] = ((i == indexes[j] ? 1 : 0) - sum) / matrix_L[j][j];
		}

		for (int j = n - 1; j >= 0; --j) {
			double sum = 0;
			for (int k = j + 1; k < n; ++k) {
				sum += matrix_U[j][k] * matrix_inverseA[k][i];
			}

			matrix_inverseA[j][i] = y[j] - sum;
		}
	}

	delete[] y;
}

double find_n1(double** matrix, int n) {
	double result = -1;
	for (int i = 0; i < n; ++i) {
		double sum = 0;
		for (int j = 0; j < n; ++j) {
			sum += fabs(matrix[i][j]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

double find_n2(double** matrix, int n) {
	double result = -1;
	for (int i = 0; i < n; ++i) {
		double sum = 0;
		for (int j = 0; j < n; ++j) {
			sum += fabs(matrix[j][i]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

double find_n3(double** matrix, int n) {
	double** tmp_matrix = nullptr;
	allocate_matrix(tmp_matrix, n);
	for (int i = 0; i < n; ++i) {
		for (int k = 0; k < n; ++k) {
			for (int j = 0; j < n; ++j) {
				tmp_matrix[i][k] += matrix[j][i] * matrix[j][k];
			}
		}
	}

	double** diagonal_matrix = get_diagonal_matrix(tmp_matrix, n);
	double result = fabs(diagonal_matrix[0][0]);
	for (int i = 1; i < n; ++i) {
		if (result < fabs(diagonal_matrix[i][i])) {
			result = fabs(diagonal_matrix[i][i]);
		}
	}

	delete_matrix(tmp_matrix, n);
	delete_matrix(diagonal_matrix, n);
	return sqrt(result);
}

double** get_diagonal_matrix(double** matrix, int n) {
	int im = 0;
	int jm = 0;
	double** new_matrix = nullptr;
	allocate_matrix(new_matrix, n);
	copy_matrix(new_matrix, matrix, n);

	double** tmp_matrix = nullptr;

	allocate_matrix(tmp_matrix, n);

	while (max_not_diagonal(new_matrix, n, im, jm) > EPS) {
		double alpha = atan(2 * new_matrix[im][jm] / (new_matrix[im][im] - new_matrix[jm][jm])) / 2;

		double** matrix_T = nullptr;
		allocate_matrix(matrix_T, n);
		for (int i = 0; i < n; ++i) {
			matrix_T[i][i] = 1;
		}

		matrix_T[im][im] = cos(alpha);
		matrix_T[jm][jm] = cos(alpha);
		matrix_T[im][jm] = -sin(alpha);
		matrix_T[jm][im] = sin(alpha);

		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < n; ++k) {
				tmp_matrix[i][k] = 0;
				for (int j = 0; j < n; ++j) {
					tmp_matrix[i][k] += matrix_T[j][i] * new_matrix[j][k];
				}
			}
		}

		delete_matrix(new_matrix, n);
		new_matrix = multiply_matrix(tmp_matrix, matrix_T, n);

		delete_matrix(matrix_T, n);
	}

	delete_matrix(tmp_matrix, n);
	return new_matrix;
}

double max_not_diagonal(double** matrix, int n, int& i, int& j) {
	double res = -1;
	for (int k = 0; k < n; ++k) {
		for (int p = 0; p < n; ++p) {
			if (k == p) continue;

			if (res < fabs(matrix[k][p])) {
				res = fabs(matrix[k][p]);
				i = k;
				j = p;
			}
		}
	}

	return res;
}

double** trans_matrix(double** matrix, int n) {
	double** new_matrix = nullptr;
	allocate_matrix(new_matrix, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			new_matrix[i][j] = matrix[j][i];
		}
	}

	return new_matrix;
}