#include "Player.hpp"
#include <cstdlib>
#include <iostream>

#include <algorithm> // sort()
#include "hmm.h"
//#include <chrono> //time
#include <float.h>

using namespace std;

//#define VERBOSE
//#define LOCAL
#define KATTIS

namespace ducks
{

#ifdef LOCAL
	static const double A_matrix_focus_percentage = 0.9;
	static const double B_matrix_focus_percentage = 0.0;

	static const double SHOOT_THRESHOLD = 0.5;
	static const int MIN_SAMPLE_SIZE_FOR_CONSIDERATION = 500;
#endif
#ifdef KATTIS
	static const double A_matrix_focus_percentage = 0.9;
	static const double B_matrix_focus_percentage = 0.0;

	static const double SHOOT_THRESHOLD = 0.6;
	static const int MIN_SAMPLE_SIZE_FOR_CONSIDERATION = 1500;
#endif

	Player::Player()
		:
		c_N(5),
		c_M(9),
		c_default_A({
		HMM::similar_stochastic_vector_focused(Player::c_N, { 0 }, A_matrix_focus_percentage),
		HMM::similar_stochastic_vector_focused(Player::c_N,{ 1 }, A_matrix_focus_percentage),
		HMM::similar_stochastic_vector_focused(Player::c_N,{ 2 }, A_matrix_focus_percentage),
		HMM::similar_stochastic_vector_focused(Player::c_N,{ 3 }, A_matrix_focus_percentage),
		HMM::similar_stochastic_vector_focused(Player::c_N,{ 4 }, A_matrix_focus_percentage) }),
		c_default_B({
		HMM::similar_stochastic_vector(Player::c_M),
		HMM::similar_stochastic_vector(Player::c_M),
		HMM::similar_stochastic_vector(Player::c_M),
		HMM::similar_stochastic_vector(Player::c_M),
		HMM::similar_stochastic_vector(Player::c_M) }),
		c_default_pi({ HMM::similar_stochastic_vector(Player::c_N) }),
		m_all_species(static_cast<int>(ESpecies::COUNT_SPECIES), HMM(Player::c_default_A, Player::c_default_B, Player::c_default_pi))
	{
		c_default_pi.print();
		c_default_A.print();
		c_default_B.print();
	}

	struct ProbabilityAction {
		ESpecies i_species;
		double p_species;
		double p_action;
		int data_set_size;
		Action action;
		vector<int> obsv_sequence;

		bool operator < (const ProbabilityAction& pa) const
		{
			return (p_species - p_action > pa.p_species - pa.p_action);
		}
	} typedef ProbabilityAction;

	void Player::train_all_species_hmm(const Deadline& p_due)
	{
		for (int i_species = 0; i_species < m_all_species.size(); ++i_species)
		{
			HMM& l_species = m_all_species[i_species];

			if (l_species.get_len_obsv_sequence() > l_species.get_len_trained_obsv_sequence()) {
				bool completed = l_species.baum_welch_algorithm(p_due);
				if (!completed)
					return;
			}
		}
	}

	/*
	Get the most probable species in species vector
	*/
	int Player::max_probable_species_for_moves(const vector<int>& obsv_sequence, double& p_max_species) const
	{
		int max_species{ static_cast<int>(SPECIES_UNKNOWN) };
		p_max_species = -DBL_MAX;

		if (obsv_sequence.size() == 0) max_species;

		for (int i = 0; i < m_all_species.size(); ++i)
		{
			const HMM& species_hmm = m_all_species[i];
			double p_species = species_hmm.log_probability_of_obsv_sequence(obsv_sequence);

			if (p_species > p_max_species) {
				p_max_species = p_species;
				max_species = i;
			}
		}
		return max_species;
	}

	vector<int> Player::get_bird_movement_sequence(const Bird& bird) const
	{
		vector<int> movements(bird.getSeqLength());
		for (int i = 0; i < bird.getSeqLength(); ++i) {
			movements[i] = bird.getObservation(i);
		}
		return movements;
	}

	vector<int> Player::get_bird_movement_sequence_no_deaths(const Bird& bird) const
	{
		vector<int> movements(bird.getSeqLength());
		int obsv;
		for (int i = 0; i < bird.getSeqLength(); ++i) {
			obsv = bird.getObservation(i);
			if (obsv == EMovement::MOVE_DEAD)
				break;
			movements[i] = obsv;
		}
		return movements;
	}

	double Player::species_trained_percentage()
	{
		int n_species = m_all_species.size();
		int n_species_trained = 0;
		for (int i_species = 0; i_species < m_all_species.size(); ++i_species)
		{
			if (m_all_species[i_species].get_len_trained_obsv_sequence() > 0) {
				n_species_trained += 1;
			}
		}
		double p_species_trained = (double)n_species_trained / (double)n_species;
		return p_species_trained;
	}

	static int round_shots = 0;
	static int round_hit_shots = 0;
	static int total_shots = 0;
	static int total_hit_shots = 0;

	Action Player::shoot(const GameState &pState, const Deadline &pDue)
	{
		int step = pState.getBird(0).getSeqLength();
		if (step == 1) {
			round_shots = 0;
			round_hit_shots = 0;
		}

		vector<ProbabilityAction> max_p_action_per_bird;

		int min_seq_len = 100 - pState.getNumBirds();
		if (step > min_seq_len && species_trained_percentage() == 1.0)
		{
			for (int i_bird = 0; i_bird < pState.getNumBirds(); ++i_bird)
			{
				const Bird& bird = pState.getBird(i_bird);

				if (bird.isAlive())
				{
					vector<int> moves = get_bird_movement_sequence(bird);
					double p_max_probable_species{ 0.0 };
					int max_probable_species = max_probable_species_for_moves(moves, p_max_probable_species);

					if (max_probable_species != ESpecies::SPECIES_UNKNOWN &&
						m_all_species[max_probable_species].get_len_trained_obsv_sequence() > MIN_SAMPLE_SIZE_FOR_CONSIDERATION &&
						max_probable_species != ESpecies::SPECIES_BLACK_STORK)
					{
						HMM& species_hmm = m_all_species[max_probable_species];

						int max_next_move{ 0 };
						double p_max_next_move{ 0.0 };
						species_hmm.max_probable_next_observation(max_next_move, p_max_next_move, moves);

						max_p_action_per_bird.push_back(
							ProbabilityAction{ static_cast<ESpecies>(max_probable_species), p_max_probable_species, p_max_next_move, species_hmm.get_len_trained_obsv_sequence(), Action{ i_bird, EMovement(max_next_move) }, moves });
					}
				}
			}
		}

		if (max_p_action_per_bird.size() > 0)
		{
			////Sort to ascending order of probability, sort function in struct
			sort(max_p_action_per_bird.begin(), max_p_action_per_bird.end());

			ProbabilityAction max_p_action = max_p_action_per_bird.back();

			train_all_species_hmm(pDue);

			if (max_p_action.p_species - max_p_action.p_action < SHOOT_THRESHOLD) {
				++round_shots;
				++total_shots;
				return max_p_action.action;
			}
			return cDontShoot;
		}

		train_all_species_hmm(pDue);
		return cDontShoot;
	}

	void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
	{
		++round_hit_shots;
		++total_hit_shots;
		std::cerr << "HIT BIRD!!!" << std::endl;
	}

	std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
	{

		cerr << "total_shots_hit:" << total_hit_shots << endl;
		cerr << "total_shots_made:" << total_shots << endl;
		cerr << "round_hit_percentage:" << round_hit_shots / (double)round_shots << endl;
		cerr << "total_hit_percentage:" << total_hit_shots / (double)total_shots << endl;
		cerr << "guess!" << endl;

		vector<ESpecies> l_guesses(pState.getNumBirds(), SPECIES_PIGEON);

		for (int i_bird = 0; i_bird < pState.getNumBirds(); ++i_bird)
		{
			const Bird& bird = pState.getBird(i_bird);
			vector<int> moves = get_bird_movement_sequence_no_deaths(bird);

			double p_max_probable_species{ 0.0 }; //TODO use this
			int max_probable_species = max_probable_species_for_moves(moves, p_max_probable_species);

			if (max_probable_species >= SPECIES_PIGEON && max_probable_species < COUNT_SPECIES) {
				l_guesses[i_bird] = static_cast<ESpecies>(max_probable_species);
			}

		}

		train_all_species_hmm(pDue);

		return l_guesses;
	}

	void Player::add_bird_moves_to_species_hmm(const GameState& p_state, const vector<ESpecies>& p_species)
	{
		for (int i_p_species = 0; i_p_species < p_species.size(); ++i_p_species)
		{
			const ESpecies& species = p_species[i_p_species];
			HMM& species_hmm = m_all_species[static_cast<int>(species)];
			if (

#ifdef LOCAL
				species_hmm.get_len_obsv_sequence() < 1000  //max my pc can handle
#endif
#ifdef KATTIS
				true
#endif
				)
			{
				const Bird& bird = p_state.getBird(i_p_species);
				vector<int> i_bird_moves = get_bird_movement_sequence_no_deaths(bird);
				species_hmm.add_to_obsv_sequence(i_bird_moves);
			}
		}
	}

	void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
	{
		add_bird_moves_to_species_hmm(pState, pSpecies);
	}


} /*namespace ducks*/
