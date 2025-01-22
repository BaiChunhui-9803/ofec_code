/********* Begin Register Information **********
{
	"name": "HGHE-DE",
	"identifier": "HGHE_DE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/***********************************************
@article{wang2023history,
  volume = {27},
  number = {6},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {1962--1975},
  year = {2023},
  author = {Junchen Wang and Changhe Li and Sanyou Zeng and Shengxiang Yang},
  title = {History-guided hill exploration for evolutionary computation}
}
***********************************************/

#ifndef OFEC_HGHE_DE_H
#define OFEC_HGHE_DE_H

#include "../../../../template/framework/hghe/continuous/hghec.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class HGHE_DE : virtual public HGHEC {
		OFEC_CONCRETE_INSTANCE(HGHE_DE)
	protected:
		enum class PtnlPredStrat { kMean, kRandom, kPareto, kObjFirst, kFreqFirst, kNumOpts };
		PtnlPredStrat m_ptnl_pred_strat;
		Real m_scaling_factor, m_crossover_rate;
		size_t m_size_subpop_exploit;
		size_t m_num_sols_explore;
		MultiPopulation<PopulationDE<>> m_subpops_exploit;
		std::vector<std::unique_ptr<Solution<>>> m_sols_explore;
		Real m_sum_potential_explore;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		int classifyPopulation(Environment *env);
		void updatePotentialExplore(Environment *env);
		int exploreHills(Environment *env);
		int exploitHills(Environment *env);

	private:
		Hill* rouletteWheelSelection();
	};
}

#endif // !OFEC_HGHE_DE_H
