#ifndef AMP_CONTINUOUS_H
#define AMP_CONTINUOUS_H

#include "../amp.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "amp_cont_pop.h"
#include "hibernated_area.h"

namespace ofec {
	template <typename TPopulation>
	class ContAMP : public AMP<PopContAMP<TPopulation>> {
	public:
		using typename AMP<PopContAMP<TPopulation>>::IndividualType;
	protected:
		Real m_conv_thresh = 0.0001;
		Real m_conv_factor = 0.005;
		HibernatedArea m_hibernated_area;
		
		void initSolutions(int Solution_num, Problem *pro, Algorithm *alg, Random *rnd) override {
			if (this->m_new_indis.size() > Solution_num) {
				this->m_new_indis.resize(Solution_num);
			}
			else {
				int originSize = this->m_new_indis.size();
				this->m_new_indis.resize(Solution_num);
				size_t number_objectives = pro->numberObjectives();
				size_t num_cons = pro->numberConstraints();
				size_t num_vars = CAST_CONOP(pro)->numberVariables();
				for (int i(originSize); i < Solution_num; ++i) {
					this->m_new_indis[i].reset(new IndividualType(number_objectives, num_cons, num_vars));
					this->m_new_indis[i]->initialize(i, pro, rnd);
					this->m_new_indis[i]->evaluate(pro, alg);
				}
			}
		}
		void updateStagnateState(Problem *pro) override;
		void removeRedundentHiber(Problem *pro) override;
		bool checkIncreaseDiv(Problem *pro) override;
		void degradeExploredAreas(Problem *pro, Algorithm *alg, Random *rnd) override;
		bool judgeOverlapping(const PopContAMP<TPopulation>& p1, const PopContAMP<TPopulation>& p2, Problem *pro) override {
			Real center_dis(pro->variableDistance(p1.at(p1.bestIdx(pro)), p2.at(p2.bestIdx(pro))));
			Real dis(0);
			if (center_dis < p1.initialRadius() && center_dis < p2.initialRadius()) {
				int c1(0), c2(0);
				for (size_t k = 0; k < p2.size(); k++) {
					dis = pro->variableDistance(p1.center(), p2[k]);
					if (dis < p1.initialRadius()) c1++;
				}

				for (size_t k = 0; k < p1.size(); k++) {
					dis = pro->variableDistance(p2.center(), p1[k]);
					if (dis < p2.initialRadius()) c2++;
				}
				return c1 > 0 && c2 > 0;
			}
			else
				return false;
		}

		void hibernatePop(PopContAMP<TPopulation> &cur_pop, Problem *pro, Algorithm *alg, Random *rnd) override {
			AMP<PopContAMP<TPopulation>>::hibernatePop(cur_pop, pro, alg, rnd);
			int idx(0);
			Real minDis(0);
			for (auto& it : cur_pop.best()) {
				m_hibernated_area.get_nearest_optimum(*it, idx, minDis, pro);
				if (minDis > CAST_CONOP(pro)->variableAccuracy()) {
					m_hibernated_area.add_optimum(*it, pro, alg, rnd);
				}
			}
		}
	
	public:
		ContAMP(size_t size_pop);
	};

	template<typename TPopulation>
	void ContAMP<TPopulation>::updateStagnateState(Problem *pro){
		for (size_t i = 0; i < this->m_pops.size(); ++i) {
			if (this->m_pops[i]->isActive()) {
                this->m_pops[i]->updateStagnant(this->m_avgRadius, pro);
			}
		}
	}

	template<typename TPopulation>
	void ContAMP<TPopulation>::removeRedundentHiber(Problem *pro) {
		int idx(0);
		Real minDis(0);
		for (auto iter = this->m_pops.begin(); iter != this->m_pops.end();) {
			if (!(*iter)->isActive()) {
				m_hibernated_area.get_nearest_optimum(*((*iter)->best(pro).front()), idx, minDis, pro);
				if (minDis < CAST_CONOP(pro)->variableAccuracy())
					iter = this->m_pops.erase(iter);
				else
					iter++;
			}
			else
				iter++;
		}
	}

	template<typename TPopulation>
	bool ContAMP<TPopulation>::checkIncreaseDiv(Problem *pro) {
		if (this->m_pops.size() == 0)
			return true;

		Real avgRadius = 0;
		int count = 0;
		for (unsigned k = 0; k < this->m_pops.size(); k++) {
			if (this->m_pops[k]->isActive() && !this->m_pops[k]->isStagnant()) {
				avgRadius += this->m_pops[k]->curRadius();
				count++;
			}
		}
		if (count > 0)	avgRadius /= count;
		else return true;
		if (avgRadius <= this->m_conv_factor * CAST_CONOP(pro)->domainArea()) {//
			return true;
		}
		else return false;
	}

	template<typename TPopulation>
	void ContAMP<TPopulation>::degradeExploredAreas(Problem *pro, Algorithm *alg, Random *rnd) {
		for (int i(0); i < this->m_pops.size(); ++i) {
			this->m_pops[i]->degradeExploredAreas(m_hibernated_area, pro, alg, rnd);
		}
	}

	template<typename TPopulation>
	ContAMP<TPopulation>::ContAMP(size_t size_pop) :
		AMP<PopContAMP<TPopulation>>(size_pop) {}
}

#endif