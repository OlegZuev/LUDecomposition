#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void allocate_matrix(double**& matrix, int& n);
void delete_matrix(double** matrix, int& n);
void copy_matrix(double** matrix_to, double** matrix_from, int n);
void print_matrix(double** matrix, int n, ostream& ostr);
double** multiply_matrix(double** matrix_a, double** matrix_b, int n);
double** subtraction_matrix(double** matrix_a, double** matrix_b, int n);

int main() {
	double** matrix_A = nullptr;
	double** matrix_PA = nullptr;
	double** matrix_inverseA = nullptr;
	double** matrix_U = nullptr;
	double** matrix_L = nullptr;
	int* indexes = nullptr;

	int n;
	ifstream fin("../input.txt");
	fin >> n;

	ostream& out = cout;

	indexes = new int[n];
	for (int i = 0; i < n; ++i) {
		indexes[i] = i;
	}

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

	print_matrix(matrix_PA, n, out);

	int permutation_count = 0;
	for (int i = 0; i < n; ++i) {
		double column_max = matrix_U[i][i];
		int k = i;
		for (int j = 0; j < n; ++j) {
			if (abs(column_max) < abs(matrix_U[j][i])) {
				column_max = matrix_U[j][i];
				k = j;
			}
		}

		permutation_count++;
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

		print_matrix(matrix_U, n, out);

		for (int j = 0; j <= i; ++j) {
			double vector_mult_result = 0;
			for (int k = 0; k < i; ++k) {
				vector_mult_result += matrix_L[i][k] * matrix_U[k][j];
			}

			matrix_L[i][j] = matrix_PA[i][j] - vector_mult_result;
		}

		print_matrix(matrix_L, n, out);
		out << "permutation count: " << permutation_count << endl;
	}

	for (int i = 0; i < n; i++) {
		out << indexes[i] << " ";
	}
	out << endl;

	double** tmp = multiply_matrix(matrix_L, matrix_U, n);
	double** error = subtraction_matrix(tmp, matrix_PA, n);
	delete_matrix(tmp, n);

	out << "L*U - P*A" << endl;
	print_matrix(error, n, out);

	double det_A = 1;
	for (int i = 0; i < n; i++) {
		det_A *= matrix_L[i][i];
	}

	det_A *= permutation_count % 2 == 0 ? 1 : -1;

	out << "Determinant=" << det_A << endl;

	print_matrix(matrix_A, n, out);

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


	//tmp = multiply_matrix(matrix_U, matrix_L, n);
	//for (int i = 0; i < n; ++i) {
	//	for (int j = 0; j < n; ++j) {
	//		tmp[i][j] /= det_A;
	//	}
	//}


	print_matrix(matrix_inverseA, n, out);
	// delete_matrix(tmp, n);

	delete_matrix(matrix_PA, n);
	delete_matrix(matrix_inverseA, n);
	delete_matrix(matrix_A, n);
	delete_matrix(matrix_U, n);
	delete_matrix(matrix_L, n);
	delete_matrix(error, n);
	delete[] indexes;
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
