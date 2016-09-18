#include "Player.hpp"
#include "Client.hpp"
#include "GameServer.hpp"

#include <iostream>
#include <fstream>

#include "hmm.h"
#include "util.h"

#include "Deadline.hpp"
#include <algorithm>

bool gVerbose = false;

void read_sequence(vector<int>& obsv_sequence)
{
	string line;
	getline(cin, line);

	stringstream ss(line);

	int value;
	while (ss >> value) {
		obsv_sequence.push_back(value);
	}
}

int main(int argc, char **argv)
{

	////////////////TEST
	/*vector<int> obsv_sequence;
	read_sequence(obsv_sequence);
	cerr << "read_done" << endl;

	obsv_sequence.insert(obsv_sequence.end(), obsv_sequence.begin(), obsv_sequence.end());

	const int c_N = 5;
	const int c_M = 9;

	Matrix<double> c_default_A({ HMM::similar_stochastic_vector_focused(c_N,{ 0 }, 0.99), HMM::similar_stochastic_vector_focused(c_N,{ 1 }, 0.99), HMM::similar_stochastic_vector_focused(c_N,{ 2 }, 0.99), HMM::similar_stochastic_vector_focused(c_N,{ 3 }, 0.99), HMM::similar_stochastic_vector_focused(c_N,{ 4 }, 0.99) });
	Matrix<double> 	c_default_B({ HMM::similar_stochastic_vector(c_M), HMM::similar_stochastic_vector(c_M), HMM::similar_stochastic_vector(c_M), HMM::similar_stochastic_vector(c_M), HMM::similar_stochastic_vector(c_M) });
	Matrix<double> 	c_default_pi({ HMM::similar_stochastic_vector(c_N) });


	HMM hmm{ c_default_A, c_default_B, c_default_pi };

	print_vector(obsv_sequence);
	hmm.add_to_obsv_sequence(obsv_sequence);
	bool status = hmm.baum_welch_algorithm(ducks::Deadline(999999));


	vector<int> sub_obsv_seq(obsv_sequence.begin() + 20, obsv_sequence.begin() + 60);
	vector<int> test_obsv_seq{ 1, 0, 4, 0, 0, 0, 4, 0, 0, 0, 1, 3, 0, 0, 3, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 7, 0, 0, 6, 7, 7, 1, 0, 7, 1, 7, 0, 6, 0, 6, 7, 6, 7, 7, 1, 6, 6, 6, 1, 0, 7, 1, 1, 6, 7, 7, 1, 7, 7, 2, 0, 0, 0, 0, 1, 0, 3, 0, 0, 3, 3, 0, 0, 0, 3, 1, 0, 0, 3, 0, 0, 3, 3, 3, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 3, };
	for (int i = 0; i < 16; ++i)
	{
		int max_obsv;
		double p;
		random_shuffle(test_obsv_seq.begin(), test_obsv_seq.end());
		hmm.log_probability_of_obsv_sequence(test_obsv_seq);
	}*/
	///////////////////

	// Parse parameters
	bool lCreateServer = false;
	std::string lLoadFilename = "ParadiseEmissions.in";

	for (int i = 1; i < argc; ++i)
	{
		std::string param(argv[i]);
		if (param == "server" || param == "s")
		{
			lCreateServer = true;
		}
		else if (param == "verbose" || param == "v")
		{
			gVerbose = true;
		}
		else if (param == "load" || param == "l")
		{
			++i;
			if (i < argc)
				lLoadFilename = argv[i];
			else
			{
				std::cerr << "Observations file must be given as an argument" << std::endl;
				exit(-1);
			}
		}
		else
		{
			std::cerr << "Unknown parameter: '" << argv[i] << "'" << std::endl;
			exit(-1);
		}
	}

	/**
	 * Start the program either as a server or a client
	 */
	if (lCreateServer)
	{
		// Create a server
		ducks::GameServer lGameServer(std::cin, std::cout);

		if (!lLoadFilename.empty())
		{
			if (gVerbose)
				std::cerr << "Loading '" << lLoadFilename << "'" << std::endl;
			std::ifstream lFile(lLoadFilename);
			lGameServer.load(lFile);
		}

		// Run the server
		lGameServer.run();
	}
	else
	{
		// Create the player
		ducks::Player lPlayer;

		// Create a client with the player
		ducks::Client lClient(lPlayer, std::cin, std::cout);

		// Run the client
		lClient.run();
	}
	return 0;
}
