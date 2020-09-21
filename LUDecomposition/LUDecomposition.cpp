#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void allocate_matrix(double**& matrix, int& n);
void delete_matrix(double** matrix, int& n);
void copy_matrix(double** matrix_to, double** matrix_from, int n);
void print_matrix(double** matrix, int n, ostream& ostr);

int main() {
	double** matrix_A = nullptr;
	double** matrix_U = nullptr;
	double** matrix_L = nullptr;

	int n;
	ifstream fin("../input.txt");
	fin >> n;

	allocate_matrix(matrix_A, n);
	allocate_matrix(matrix_U, n);
	allocate_matrix(matrix_L, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> matrix_A[i][j];
		}
	}

	fin.close();
	copy_matrix(matrix_U, matrix_A, n);

	print_matrix(matrix_A, n, cout);

	for (int i = 0; i < n; ++i) {
		double column_max = matrix_U[i][i];
		int k = i;
		for (int j = 0; j < n; ++j) {
			if (abs(column_max) < abs(matrix_U[j][i])) {
				column_max = matrix_U[j][i];
				k = j;
			}
		}

		for (int j = 0; j < n; ++j) {
			swap(matrix_U[i][j], matrix_U[k][j]);
			swap(matrix_A[i][j], matrix_A[k][j]);
		}

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
		
		print_matrix(matrix_U, n, cout);

		for (int j = 0; j <= i; ++j) {
			double vector_mult_result = 0;
			for (int k = 0; k < i; ++k) {
				vector_mult_result += matrix_L[i][k] * matrix_U[k][j];
			}

			matrix_L[i][j] = matrix_A[i][j] - vector_mult_result;
		}

		print_matrix(matrix_L, n, cout);
	}

	delete_matrix(matrix_A, n);
	delete_matrix(matrix_U, n);
	delete_matrix(matrix_L, n);
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
