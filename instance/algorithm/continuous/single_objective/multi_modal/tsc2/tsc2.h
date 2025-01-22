/********* Begin Register Information **********
{
    "description": "Topological species sonservation algorithm",
    "identifier": "TSC2",
    "name": "TSC2",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/***********************************************
@article{stoean2010multimodal,
  volume = {14},
  number = {6},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {842--864},
  year = {2010},
  author = {Catalin Stoean and Mike Preuss and Ruxandra Stoean and Dumitru Dumitrescu},
  title = {Multimodal optimization by means of a topological species conservation algorithm}
}
***********************************************/

#ifndef OFEC_TSC2_H
#define OFEC_TSC2_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../core/algorithm/population.h"

namespace ofec {
	class IndTSC2 : public Solution<> {
	public:
		bool marked = false, marked_as_seeds = false, selected_for_crossover = false, found = false;
		int former_index = 0, former_id = 0, seed_id = 0, old_seed_id = 0;
		std::vector<int> ids;

		IndTSC2(size_t number_objectives, size_t num_cons, size_t num_vars) : 
			Solution<>(number_objectives, num_cons, num_vars) {}
	};

	class TSC2 : virtual  public Algorithm {
		OFEC_CONCRETE_INSTANCE(TSC2)
	protected:
		size_t m_pop_size, m_num_gradations, m_tournament_size;
		Real m_percent_seeds;
		Real m_pc, m_pm, m_ms;
		Population<IndTSC2> m_pop, m_seeds;
		size_t m_iteration;
		std::vector<Real> m_gradations;
		Real m_positive;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		void determineSeeds(Environment *env);
		void tournamentSelection(Environment *env);
		void recombination(Environment *env);
		void normalMutation(Environment *env);
		void insertSeedsWithoutDuplicates(Environment *env);
		void integrateSolutionsToExistingSeeds(Environment *env);
		size_t countFreeSolutions();
		void integrateFreeSolutions(Environment *env);
		bool containsSeed(const IndTSC2 &ind, int seed_id);
		bool hillValley(const IndTSC2 &ind1, const IndTSC2 &ind2, Environment *env);
		void assignIDs(Environment *env);
	};
}

#endif // ! OFEC_TSC2_H
