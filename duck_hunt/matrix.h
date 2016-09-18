#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

using namespace std;

template <typename T>
class Matrix {
public:
	vector<vector<T>> matrix;

	Matrix() : matrix(vector<vector<T>>()) {
	}

	Matrix(int rows, int cols) : matrix(vector<vector<T>>(rows, vector<T>(cols))) {
	}

	Matrix(vector<vector<T>> p_matrix) : matrix(p_matrix)
	{

	}

	static Matrix<T> multiply(Matrix<T> A, Matrix<T> B) {
		T sum = 0;
		int c, d, k;
		Matrix<T> C = Matrix<T>(A.matrix.size(), B.matrix[1].size());
		int m = A.matrix.size();
		int n = A.matrix[0].size();
		int p = B.matrix.size();
		int q = B.matrix[0].size();
		if (n == p) {
			for (c = 0; c < m; c++) {
				for (d = 0; d < q; d++) {
					for (k = 0; k < p; k++) {
						sum = sum + A.matrix[c][k] * B.matrix[k][d];
					}
					C.matrix[c][d] = sum;
					sum = 0;
				}
			}
		}
		else {
			cerr << "Mismatched matrix dimensions" << endl;
		}
		return C;
	}

	void mul_to_the_right(Matrix<T> a) {
		this->matrix = multiply(this->matrix, a);
	}

	void mul_to_the_left(Matrix<T> a) {
		this->matrix = multiply(a, this->matrix);
	}

	static Matrix<T> transpose(Matrix<T> A) {
		return 0;
	}

	void print() const {
		for (unsigned int i = 0; i < this->matrix.size(); i++) {
			for (unsigned int j = 0; j < this->matrix[i].size(); j++)
				cerr << this->matrix[i][j] << " ";
			cerr << endl;
		}
		cerr << endl;
	}

	void print_inline() const {
		cerr << this->matrix.size() << " " << this->matrix[0].size() << " ";
		for (int i = 0; i < this->matrix.size(); i++) {
			for (int j = 0; j < this->matrix[i].size(); j++)
				cerr << this->matrix[i][j] << " ";
		}
		cerr << endl;
	}
};

#endif //MATRIX_H