#ifndef AMP_POP_CONTINUOUS_H
#define AMP_POP_CONTINUOUS_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/global.h"
#include "hibernated_area.h"
#include <vector>
#include "../amp_pop.h"

namespace ofec {
	template <typename TPopulation>
	class PopContAMP : public PopAMP<TPopulation> {
	public:
		using IndividualType = typename TPopulation::IndividualType;

	protected:
		Real m_converg_threshold = 0.0001;
		Real m_converg_factor = 0.005;
		IndividualType m_center;

	public:
		PopContAMP(size_t pop_size, Problem *pro) : 
			PopAMP<TPopulation>(pop_size, pro), 
			m_center(pro->numberObjectives(), 
				pro->numberConstraints(),
				CAST_CONOP(pro)->numberVariables()) {}

		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
			int rf = PopAMP<TPopulation>::evolve(pro, alg, rnd);
			this->updateBest(pro);
			if (rf)return rf;
			std::vector<IndividualType> new_inds;
			for (auto& best : this->m_best) {
				new_inds.push_back(*best);
				if (this->m_stagnant) {
					new_inds.back().cauchyMove(pro, rnd);
				}
				else {
					new_inds.back().brwonianMove(pro, rnd, this->m_cur_radius);
				}
			}
			for (auto& new_ind : new_inds) {
				rf = new_ind.evaluate(pro, alg);
				this->updateBest(new_ind, pro);
				if (rf )return rf;
			}
			return rf;
		}

		void setParameters(Real convergThreshold, Real convFactor) {
			m_converg_threshold = convergThreshold;
			m_converg_factor = convFactor;
		}

		virtual void degradeExploredAreas(HibernatedArea& h, Problem *pro, Algorithm *alg, Random *rnd) {
			for (auto& it : this->m_best) {
				h.derateFitness(*it, pro, alg, rnd);
			}
			for (auto& it : this->m_individuals) {
				h.derateFitness(*it, pro, alg, rnd);
				// for pso
				//h.derateFitness(it->pbest());
			}
		}

		void updateStagnant(double avgRadius, Problem *pro) {
			if (this->m_count_radius_decre >= 1. * this->size() && 
				this->m_cur_radius >= avgRadius && 
				this->m_cur_radius > m_converg_factor * CAST_CONOP(pro)->domainArea() || 
				this->m_count_radius_decre >= 10 * this->size())
				this->m_stagnant = true;
			else
				this->m_stagnant = false;
		}

		bool judgeConverged() override {
			return this->m_cur_radius <= m_converg_threshold;
		}

		void udpateCurRadius(Problem *pro) override {
			calculateCenter();
			this->m_cur_radius = this->avgDisTo(m_center, pro);
		}

		const IndividualType& center() const {
			return m_center;
		}

		void calculateCenter() {
			if (this->m_individuals.empty())return;
			m_center = this->m_individuals.front()->phenotype();
			size_t var_size(m_center.variable().size());
			for (decltype(this->m_individuals.size()) i(1); i < this->m_individuals.size(); ++i) {
				for (int j(0); j < var_size; ++j) {
					m_center.variable()[j] += this->m_individuals[i]->variable()[j];
				}
			}
			for (int j(0); j < var_size; ++j) {
				m_center.variable()[j] /= this->m_individuals.size();
			}
		}

		int bestIdx(Problem *pro) const {
			int l = 0;
			for (int i = 0; i < this->m_individuals.size(); ++i) {
				if (this->m_individuals[i]->dominate(*(this->m_individuals[l]), pro)) {
					l = i;
				}
			}
			return l;
		}
	};
}

#endif