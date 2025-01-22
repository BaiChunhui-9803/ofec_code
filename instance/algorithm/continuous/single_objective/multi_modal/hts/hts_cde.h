/********* Begin Register Information **********
{
	"identifier": "HTS_CDE",
	"name": "HTS-CDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
@article{li2015history,
  volume = {19},
  number = {1},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {136--150},
  year = {2015},
  author = {Lingxi Li and Ke Tang},
  title = {History-based topological speciation for multimodal optimization}
}
******************************************************************************/

#ifndef OFEC_HTS_CDE_H
#define OFEC_HTS_CDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../crowding_de/crowding_de.h"

namespace ofec {
	class SC_CDE : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(SC_CDE)
	protected:
		size_t m_pop_size;
		std::unique_ptr<PopulationDE<>> m_pop;
		std::list<Solution<>> m_seeds;
		std::list<std::shared_ptr<const Solution<>>> m_history_points;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		void speciation(Environment *env);
		void conserveSeeds(Environment *env);
		virtual bool sameSpecies(const Solution<> &a, const Solution<> &b, 
			const std::list<std::shared_ptr<const Solution<>>> &H, Environment *env
		) = 0;
	};

	class HTS_CDE : public SC_CDE {
		OFEC_CONCRETE_INSTANCE(HTS_CDE)
	protected:
		void addInputParameters() {}
		bool sameSpecies(const Solution<> &a, const Solution<> &b, 
			const std::list<std::shared_ptr<const Solution<>>> &H, Environment *env
		) override;

	public:
		static bool test(const Solution<> &s1, const Solution<> &s2, const std::list<std::shared_ptr<const Solution<>>> &H, Environment *env);
		static std::vector<std::pair<Real, Real>> cubeRegion(const Solution<> &s1, const Solution<> &s2);
		static bool isPointInOpenCubeRegion(const Solution<> &p, const std::vector<std::pair<Real, Real>> &R);
		static std::shared_ptr<const Solution<>> mediumPoint(const Solution<> &s1, const Solution<> &s2, const std::list<std::shared_ptr<const Solution<>>> &H);
	};
}

#endif // ! OFEC_HTS_CDE_H
