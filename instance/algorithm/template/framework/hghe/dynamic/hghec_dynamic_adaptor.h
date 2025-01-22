#ifndef OFEC_HGHEC_DYNAMIC_ADAPTOR_H
#define OFEC_HGHEC_DYNAMIC_ADAPTOR_H

#include "../continuous/hghec_adaptor.h"
#include "../../../../../../core/problem/encoding.h"
#include "../../../../../../utility/kd-tree/kdtree_space.h"

namespace ofec {
	class AdaptorDynamicHGHEC : public AdaptorHGHEC {
	public:
		AdaptorDynamicHGHEC(Problem *pro, Random *rnd, Algorithm *alg, size_t init_num);
		virtual void updateHills() override;
		virtual void archiveSolution(const SolutionBase& sol, TaskHGHE task) override;
		virtual void clearPreHis();
		void calculateExplorationPotential() override;
		std::vector<const SolutionBase*> getBestSolsEachHill();
	protected:
		virtual void getHisSolsInHill(std::vector<const SolutionBase*>& his_sols, int id_hill) override;
		virtual void getHisExploreSolsInHill(std::vector<const SolutionBase*>& his_sols, int id_hill);
		virtual void assignSubspaceFitness() override;
	};

}


#endif // !OFEC_HGHEC_DYNAMIC_ADAPTOR_H