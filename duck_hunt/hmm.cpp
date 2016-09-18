#include "hmm.h"
#include <math.h>
#include <sstream>
#include <algorithm> //sort()

#include <functional> // greater()
#include <numeric> //prod()
#include <functional> //prod()
#include <float.h> //DBL_MAX
#include <iterator> //back_inserter
#include "Deadline.hpp"

HMM::HMM(const Matrix<double> p_A, const Matrix<double> p_B, const Matrix<double> p_pi)
	:
	c_starting_A(p_A),
	c_starting_B(p_B),
	c_starting_pi(p_pi),
	m_A(c_starting_A),
	m_B(c_starting_B),
	m_pi(c_starting_pi),
	m_scale(0),
	m_alpha(0, 0),
	m_beta(0, 0),
	m_digamma(0, { 0, 0 }),
	m_gamma(0, 0),
	m_trained_obsv_seq_len(0),
	m_obsv_sequence(0)
{
}

void HMM::add_to_obsv_sequence(vector<int> new_obsv_sequence)
{
	std::move(new_obsv_sequence.begin(), new_obsv_sequence.end(), std::back_inserter(m_obsv_sequence));
}

bool HMM::baum_welch_algorithm(const ducks::Deadline& p_due)
{
	const vector<int> l_obsv_sequence = m_obsv_sequence;

	const unsigned int T = m_obsv_sequence.size();
	if (T == 0)	return true;
	cerr << "dataset size:" << T << endl;

	int l_N = c_starting_A.matrix.size();
	int l_M = c_starting_B.matrix[0].size();

	Matrix<double> l_pi = c_starting_pi;
	Matrix<double> l_A = c_starting_A;
	Matrix<double> l_B = c_starting_B;

	vector<double> l_scale; //c
	Matrix<double> l_alpha;
	Matrix<double> l_beta;
	vector<Matrix<double>> l_digamma;
	Matrix<double> l_gamma;

	int64_t t_buffer = T / 5.0;
	double error_buffer = T / 100.0;

	bool optimized = false;

	double p_obsv_sequence{ -DBL_MAX };
	double new_p_obsv_sequence{ -DBL_MAX };

	while (!optimized)
	{
		fill_alpha_matrix_scaled(l_pi, l_A, l_B, l_obsv_sequence, l_alpha, l_scale); if (p_due.remainingMs() < t_buffer) { return false; }
		fill_beta_matrix_scaled(l_A, l_B, l_obsv_sequence, l_scale, l_beta); if (p_due.remainingMs() < t_buffer) { return false; }
		fill_gamma_and_digamma_matrices(l_obsv_sequence, l_A, l_B, l_alpha, l_beta, l_digamma, l_gamma); if (p_due.remainingMs() < t_buffer) { return false; }

		recalculate_matrix_pi_entries(l_gamma, l_pi); if (p_due.remainingMs() < t_buffer) { return false; }
		recalculate_matrix_A_entries(l_digamma, l_gamma, l_A); if (p_due.remainingMs() < t_buffer) { return false; }
		recalculate_matrix_B_entries(l_obsv_sequence, l_gamma, l_B); if (p_due.remainingMs() < t_buffer) { return false; }

		new_p_obsv_sequence = log_probability(l_scale); if (p_due.remainingMs() < t_buffer) { return false; }

		//cerr << "oldp:" << p_obsv_sequence << " newp:" << new_p_obsv_sequence << " size:" << T << endl;
		if (new_p_obsv_sequence - p_obsv_sequence < error_buffer) {
			m_pi = l_pi;
			m_A = l_A;
			m_B = l_B;

			m_scale = l_scale;
			m_alpha = l_alpha;
			m_beta = l_beta;
			m_digamma = l_digamma;
			m_gamma = l_gamma;

			m_trained_obsv_seq_len = l_obsv_sequence.size();

			//if (m_trained_obsv_seq_len >= 198) {
			//	cerr << "defaults:" << endl;
			//	c_starting_pi.print();
			//	c_starting_A.print();
			//	c_starting_B.print();
			//	cerr << "sec_len:" << m_trained_obsv_seq_len << endl;
			//	m_pi.print();
			//	m_A.print();
			//	m_B.print();
			//	//cerr << "gamma:" << endl;
			//	//m_gamma.print();
			//}
			optimized = true;
		}
		else {
			p_obsv_sequence = new_p_obsv_sequence;
		}
	}
	return true;
}


void HMM::fill_alpha_matrix_scaled(const Matrix<double>& p_pi, const Matrix<double>& p_A, const Matrix<double>& p_B, const vector<int>& obsv_sequence, Matrix<double>& alpha, vector<double>& scale) const
{
	const int T = obsv_sequence.size();
	const int l_N = p_A.matrix.size();

	alpha = Matrix<double>(l_N, T);
	scale = vector<double>(T);

	//compute a0(i)
	double& scale_0 = scale[0];
	scale_0 = 0.0;
	for (int i = 0; i < l_N; ++i) {
		double& alpha_current = alpha.matrix[i][0];
		alpha_current = p_B.matrix[i][obsv_sequence[0]] * p_pi.matrix[0][i];

		scale_0 += alpha_current;
	}
	//scale the a0(i)
	scale_0 = 1 / scale_0;
	for (int i = 0; i < l_N; ++i) {
		alpha.matrix[i][0] = scale_0 * alpha.matrix[i][0];
	}

	//compute at(i)
	for (int t = 1; t < T; ++t) { //go through each time step
		const int& obsv = obsv_sequence[t];
		double& scale_t = scale[t];
		scale_t = 0.0;

		for (int i = 0; i < l_N; ++i) {  //go through every state i for the time step
			double& alpha_current = alpha.matrix[i][t];
			for (int j = 0; j < l_N; ++j) { //go through each possible previous state j
				alpha_current += alpha.matrix[j][t - 1] * p_A.matrix[j][i];
			}

			alpha_current = alpha_current * p_B.matrix[i][obsv];

			scale_t += alpha_current;
		}
		////scale at(i)
		scale_t = scale_t > 0 ? (1.0 / scale_t) : DBL_MAX;
		for (int i = 0; i < l_N; ++i) {
			double& alpha_current = alpha.matrix[i][t];
			alpha_current = scale_t * alpha_current;
		}
	}
}

void HMM::fill_beta_matrix_scaled(const Matrix<double>& p_A, const Matrix<double>& p_B, const vector<int>& obsv_sequence, const vector<double>& scale_vector, Matrix<double>& beta) const
{
	const int T = obsv_sequence.size();
	const int l_N = p_A.matrix.size();

	beta = Matrix<double>(l_N, T);


	for (int i = 0; i < l_N; ++i) {
		beta.matrix[i][T - 1] = scale_vector[T - 1];
	}

	for (int t = T - 2; t >= 0; --t)
	{
		const double& scale_t = scale_vector[t];
		const int& obsv_next = obsv_sequence[t + 1];
		for (int i = 0; i < l_N; ++i) //go through every state i for the time step
		{
			double& beta_current = beta.matrix[i][t];
			beta_current = 0.0;

			const vector<double>& p_A_i = p_A.matrix[i];

			for (int j = 0; j < l_N; ++j) //go through each possible following state j
			{
				const double& probability_of_transition = p_A_i[j];
				const double& probability_of_obsv_state_next = p_B.matrix[j][obsv_next];
				const double& beta_next = beta.matrix[j][t + 1];

				beta_current += probability_of_transition * probability_of_obsv_state_next * beta_next;
			}
			// scale bt(i) with same scale factor as at(i)
			beta_current = scale_t * beta_current;
		}
	}
}

void HMM::fill_gamma_and_digamma_matrices(
	const vector<int>& obsv_sequence,
	const Matrix<double>& p_A,
	const Matrix<double>& p_B,
	const Matrix<double>& alpha,
	const Matrix<double>& beta,
	vector<Matrix<double>>& digamma,
	Matrix<double>& gamma) const
{
	const int T = obsv_sequence.size();
	const int l_N = p_A.matrix.size();

	digamma = vector<Matrix<double>>(T - 1, Matrix<double>(l_N, l_N));
	gamma = Matrix<double>(l_N, T);

	for (int t = 0; t < T - 1; ++t)
	{
		vector<vector<double>>& digamma_t = digamma[t].matrix;
		const int& obsv_next = obsv_sequence[t + 1]; // get O(t+1)
		for (int i = 0; i < l_N; ++i)
		{
			const double& alpha_i_t = alpha.matrix[i][t];
			double& gamma_t_i = gamma.matrix[i][t];
			gamma_t_i = 0.0;

			vector<double>& digamma_t_i = digamma_t[i];
			const vector<double>& l_A_i = p_A.matrix[i];

			for (int j = 0; j < l_N; ++j) {
				double& digamma_t_i_j = digamma_t_i[j];

				digamma_t_i_j =
					(alpha_i_t * l_A_i[j] * p_B.matrix[j][obsv_next] * beta.matrix[j][t + 1]);

				gamma_t_i += digamma_t_i_j;
			}
		}
	}

	//Special case for T - 1
	for (int i = 0; i < l_N; ++i) {
		gamma.matrix[i][T - 1] =
			alpha.matrix[i][T - 1];
	}
}


void HMM::recalculate_matrix_pi_entries(const Matrix<double>& p_gamma, Matrix<double>& p_pi) const
{
	const int l_N = p_pi.matrix[0].size();

	for (int i = 0; i < l_N; ++i) {
		p_pi.matrix[0][i] = p_gamma.matrix[i][0];
	}
}

void HMM::recalculate_matrix_A_entries(const vector<Matrix<double>>& p_digamma, const Matrix<double>& p_gamma, Matrix<double>& p_A) const
{
	double numerator{ 0.0 };
	double denominator{ 0.0 };

	const int l_N = p_A.matrix.size();

	for (int i = 0; i < l_N; ++i) {
		const vector<double>& gamma_i = p_gamma.matrix[i];
		for (int j = 0; j < l_N; ++j) {
			for (int t = 0; t < p_digamma.size() /*=T-1*/; ++t) {
				numerator += p_digamma[t].matrix[i][j];
				denominator += gamma_i[t];
			}
			p_A.matrix[i][j] = numerator / denominator;
			numerator = 0.0;
			denominator = 0.0;
		}
	}
}


void HMM::recalculate_matrix_B_entries(const vector<int>& obsv_sequence, const Matrix<double>& p_gamma, Matrix<double>& p_B) const
{
	double numerator{ 0.0 };
	double denominator{ 0.0 };

	const int l_N = p_B.matrix.size();
	const int l_M = p_B.matrix[0].size();

	for (int j = 0; j < l_N; ++j)
	{
		const vector<double>& gamma_j = p_gamma.matrix[j];
		for (int k = 0; k < l_M; ++k)
		{
			for (int t = 0; t < obsv_sequence.size(); ++t)
			{
				const int& obsv = obsv_sequence[t];
				const double& value = gamma_j[t];
				if (obsv == k) {
					numerator += value;
				}
				denominator += value;
			}
			p_B.matrix[j][k] = numerator / denominator;
			numerator = 0.0;
			denominator = 0.0;
		}
	}
}

void HMM::max_probable_next_observation(int& max_obsv, double& p_max_obsv, const vector<int>& p_obsv_seq) const
{
	max_obsv = -1;
	p_max_obsv = -DBL_MAX;

	vector<int> obsv_seq = p_obsv_seq;
	const int T = obsv_seq.size();
	const int l_M = m_A.matrix.size();

	obsv_seq.push_back(0.0);
	int& obsv_seq_T = obsv_seq[T];

	double p_obsv;
	for (int obsv = 0; obsv < l_M; ++obsv)
	{
		obsv_seq_T = obsv;
		p_obsv = log_probability_of_obsv_sequence(obsv_seq);
		if (p_obsv > p_max_obsv) {
			p_max_obsv = p_obsv;
			max_obsv = obsv;
		}
	}
}

void HMM::max_probable_current_state(int& max_state, double& p_max_state, const vector<int> p_obsv_seq) const
{

}


double HMM::log_probability_of_obsv_sequence(const vector<int>& obsv_sequence) const
{
	Matrix<double> l_alpha;
	vector<double> l_scale;

	fill_alpha_matrix_scaled(m_pi, m_A, m_B, obsv_sequence, l_alpha, l_scale);

	double p = log_probability(l_scale);
	return p;
}

//Mark
double HMM::log_probability(const vector<double>& c_vector)
{
	double p = 0.0;
	for (int t = 0; t < c_vector.size(); ++t) {
		p += log(c_vector[t]);
	}

	return -p;
}

vector<double> HMM::similar_stochastic_vector_focused(int N, const vector<int>& focused_indices, const double focus_percentage)
{
	vector<double> base = similar_stochastic_vector(N);
	vector<double> focused_result(N);

	double focus_base_sum = 0.0;
	for (int i = 0; i < focused_indices.size(); ++i) {
		focus_base_sum += base[focused_indices[i]];
	}

	double sum_to_spare = (1.0 - focus_base_sum) * focus_percentage;

	double focused_sum = 0.0;
	for (int i = 0; i < N - 1; ++i)
	{
		if (find(focused_indices.begin(), focused_indices.end(), i) != focused_indices.end()) {
			focused_result[i] = base[i] + (base[i] / focus_base_sum) * sum_to_spare;
		}
		else {
			focused_result[i] = base[i] - (base[i] * focus_percentage);
		}
		focused_sum += focused_result[i];
	}
	focused_result[N - 1] = 1.0 - focused_sum;

	return focused_result;
}

vector<double> HMM::similar_stochastic_vector(int N)
{
	if (N == 0)
		return vector<double>(0);

	vector<double> arr(N);

	double sum = 1.0;

	double base_value = 1.0 / (double)N;
	double diff = base_value / 33.3333333;

	double& a_0 = arr[0];
	a_0 = (((2 * sum) / (double)N) - (N - 1)*diff) / 2.0;

	double a_N = 1.0 - a_0;
	for (int i = 1; i < N - 1; ++i)
	{
		arr[i] = a_0 + i * diff;
		a_N -= arr[i];
	}
	arr[N - 1] = a_N;

	random_shuffle(arr.begin(), arr.end());

	return arr;
}

vector<double> HMM::random_stochastic_vector(int n, vector<int> prioritized_indices)
{
	bool unique = false;


	//generate an vector of unique random numbers
	vector<int> arr(n);
	while (!unique)
	{
		for (int i = 0; i < n; ++i) {
			arr[i] = rand();
		}
		sort(arr.begin(), arr.end(), greater<int>());
		if (!sorted_array_contains_equal(arr))
			unique = true;
	}

	double sum = 0;
	for (int i = 0; i < arr.size(); ++i) {
		sum += arr[i];
	}

	//calculate percentages
	vector<double> s_vector(n, 0.0);
	double last = 1.0; //for calibration
	for (int i = 0; i < n - 1; ++i) {
		double val = arr[i] / sum;
		s_vector[i] = val;
		last -= val;
	}
	s_vector[n - 1] = last;


	//prioritize indices
	vector<double> prioritized_s_vector(n, 0.0);
	for (int k = 0; k < prioritized_indices.size(); ++k) {
		prioritized_s_vector[prioritized_indices[k]] = s_vector[k];
	}
	for (int i = 0, k = prioritized_indices.size(); i < n; ++i) {
		if (prioritized_s_vector[i] == 0.0) {
			prioritized_s_vector[i] = s_vector[k];
			++k;
		}
	}
	return prioritized_s_vector;
}

bool HMM::sorted_array_contains_equal(vector<int> array)
{
	for (int i = 0; i < array.size() - 1; i++)
	{
		if (array[i] == array[i + 1])
			return true;
	}
	return false;
}
