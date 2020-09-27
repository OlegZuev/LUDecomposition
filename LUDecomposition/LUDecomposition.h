#pragma once
#include <ostream>

void allocate_matrix(double**& matrix, int& n);

void delete_matrix(double** matrix, int& n);

void copy_matrix(double** matrix_to, double** matrix_from, int n);

void print_matrix(double** matrix, int n, std::ostream& ostr);

void print_array(double* arr, int n, std::ostream& ostr);

double** multiply_matrix(double** matrix_a, double** matrix_b, int n);

double** subtraction_matrix(double** matrix_a, double** matrix_b, int n);

void find_LU(double** matrix_PA, double** matrix_U, double** matrix_L, int* indexes, int n, std::ostream& out, int& permutation_count);

void find_inverse_matrix(double** matrix_inverseA, double** matrix_U, double** matrix_L, int n, int* indexes);

void find_x(double* b, double** matrix_U, double** matrix_L, int n, int* indexes);

double find_n1(double** matrix, int n);

double find_n2(double** matrix, int n);

double find_n3(double** matrix, int n);

double max_not_diagonal(double** matrix, int n, int& i, int& j);

double** get_diagonal_matrix(double** matrix, int n);