#ifndef HMM_H
#define HMM_H

#include "matrix.h"
#include <iostream>

#include "Deadline.hpp" //Deadline

using namespace std;

class HMM {
public:

	HMM(const Matrix<double> p_A, const Matrix<double> p_B, const Matrix<double> p_pi);

	int get_len_obsv_sequence() const
	{
		return m_obsv_sequence.size();
	}
	int get_len_trained_obsv_sequence() const
	{
		return m_trained_obsv_seq_len;
	}
	vector<int> get_obsv_sequence() const
	{
		return m_obsv_sequence;
	}
	void add_to_obsv_sequence(vector<int> new_obsv_sequence);

	bool baum_welch_algorithm(const ducks::Deadline& p_due);

	void max_probable_next_observation(int& max_obsv, double& p_max_obsv, const vector<int>& p_obsv_seq) const;
	void max_probable_current_state(int& max_state, double& p_max_state, const vector<int> p_obsv_seq) const;
	double log_probability_of_obsv_sequence(const vector<int>& obsv_sequence) const;

	static vector<double> similar_stochastic_vector(int N);
	static vector<double> similar_stochastic_vector_focused(int N, const vector<int>& focus_indices, double focus_percentage);
	static vector<double> random_stochastic_vector(int n, vector<int> prioritized_indices = {});
private:

	const Matrix<double> c_starting_A;
	const Matrix<double> c_starting_B;
	const Matrix<double> c_starting_pi;

	Matrix<double> m_A;
	Matrix<double> m_B;
	Matrix<double> m_pi;

	vector<double> m_scale; //c
	Matrix<double> m_alpha;
	Matrix<double> m_beta;
	vector<Matrix<double>> m_digamma;
	Matrix<double> m_gamma;

	int m_trained_obsv_seq_len;
	vector<int> m_obsv_sequence;

	void fill_alpha_matrix_scaled(const Matrix<double>& p_pi, const Matrix<double>& p_A, const Matrix<double>& p_B, const vector<int>& obsv_sequence, Matrix<double>& alpha, vector<double>& scale) const;
	void fill_beta_matrix_scaled(const Matrix<double>& p_A, const Matrix<double>& p_B, const vector<int>& obsv_sequence, const vector<double>& scale_vector, Matrix<double>& beta) const;

	void fill_gamma_and_digamma_matrices(const vector<int>& obsv_sequence, const Matrix<double>& p_A, const Matrix<double>& p_B, const Matrix<double>& alpha, const Matrix<double>& beta, vector<Matrix<double>>& digamma, Matrix<double>& gamma) const;
	void recalculate_matrix_pi_entries(const Matrix<double>& p_gamma, Matrix<double>& p_pi) const;
	void recalculate_matrix_A_entries(const vector<Matrix<double>>& p_digamma, const Matrix<double>& p_gamma, Matrix<double>& p_A) const;
	void recalculate_matrix_B_entries(const vector<int>& obsv_sequence, const Matrix<double>& p_gamma, Matrix<double>& p_B) const;

	//Util
	static double log_probability(const vector<double>& c_vector);
	static bool sorted_array_contains_equal(vector<int> array);
};
#endif //HMM_H