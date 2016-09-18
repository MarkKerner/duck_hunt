#ifndef _DUCKS_PLAYER_HPP_
#define _DUCKS_PLAYER_HPP_

#include "Deadline.hpp"
#include "GameState.hpp"
#include "Action.hpp"
#include <vector>
#include "hmm.h"

namespace ducks
{

class Player
{
public:

	const int c_N;
	const int c_M;

	const Matrix<double> c_default_A;
	const Matrix<double> c_default_B;
	const Matrix<double> c_default_pi;

	vector<HMM> m_all_species;
    /**
     * Constructor
     * There is no data in the beginning, so not much should be done here.
     */
    Player();

	void train_all_species_hmm(const Deadline& p_due);
	int max_probable_species_for_moves(const vector<int>& obsv_sequence, double& p_max_species) const;
	vector<int> get_bird_movement_sequence(const Bird& bird) const;
	vector<int> get_bird_movement_sequence_no_deaths(const Bird& bird) const;

	double species_trained_percentage();
	void add_bird_moves_to_species_hmm(const GameState& p_state, const vector<ESpecies>& p_species);
	/**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each birds contains all past actions.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    Action shoot(const GameState &pState, const Deadline &pDue);

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    std::vector<ESpecies> guess(const GameState &pState, const Deadline &pDue);

	/**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    void hit(const GameState &pState, int pBird, const Deadline &pDue);

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    void reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue);
};

} /*namespace ducks*/

#endif
