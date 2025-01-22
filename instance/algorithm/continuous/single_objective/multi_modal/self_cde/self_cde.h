/********* Begin Register Information **********
[
	{ "identifier": "SelfCCDE", "name": "self-CCDE", "tags": [ "continuous", "single-objective" ] },
	{ "identifier": "SelfCSDE", "name": "self-CSDE", "tags": [ "continuous", "single-objective" ] }
]
*********** End Register Information **********/

/***********************************************
 @article{gao2014cluster,
  volume = {44},
  number = {8},
  journal = {IEEE Transactions on Cybernetics},
  pages = {1314--1327},
  year = {2014},
  author = {Weifeng Gao and Gary G. Yen  and Sanyang Liu},
  title = {A cluster-based differential evolution with self-adaptive strategy for multimodal optimization}
}
***********************************************/

#ifndef OFEC_SELF_CDE_H
#define OFEC_SELF_CDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class SelfCDE : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(SelfCDE)
	protected:
		size_t m_pop_size, m_cluster_size;
		PopulationDE<> m_pop;
		std::vector<std::vector<size_t>> m_clusters;
		bool m_without_crowding;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		virtual void clusteringPartition(Environment *env) = 0;
		size_t nearestInd(const Solution<> &s, Environment *env);
	};

	class SelfCCDE : public SelfCDE {
		OFEC_CONCRETE_INSTANCE(SelfCCDE)
	protected:
		void addInputParameters() {}
		void clusteringPartition(Environment *env) override;
	};

	class SelfCSDE : public SelfCDE {
		OFEC_CONCRETE_INSTANCE(SelfCSDE)
	protected:
		void addInputParameters() {}
		void clusteringPartition(Environment *env) override;
	};
}

#endif // ! OFEC_SELF_CDE_H
