#ifndef AMP_POP_H
#define AMP_POP_H

#include "../../../../../core/algorithm/multi_population.h"
#include "../../../../../core/global.h"
#include <algorithm> 
#include <vector>

namespace ofec {
	template<typename>  class AMP;

	template <typename TPopulation>
	class PopAMP : public TPopulation {
	public:
		template<typename> friend class AMP;
		using IndividualType = typename TPopulation::IndividualType;

	protected:
		bool m_stagnant = false;
		Real m_cur_radius = 0;
		Real m_initial_radius = 0;
		int m_count_radius_decre = 0;

	//	int mc_max_idxs, mc_min_idxs;
	public:
		PopAMP(size_t pop_size, Problem *pro) : TPopulation(pop_size, pro) {}
		void setInitialRadius(Real r) {
			m_initial_radius = r;
		}
		Real initialRadius() const{
			return m_initial_radius;
		}
		Real curRadius() const {
			return m_cur_radius;
		}
		bool isStagnant() const {
			return m_stagnant;;
		}
		// r1 domiante r2
		Real updateRadius(Real  r1, Real r2) {
			if (r1 > r2) {
				r1 -= (r1 - r2)*r2 / r1;
			}
			else {
				r1 += (r2 - r1)*r1 / r2;
			}

			return r1;
		}

		Real avgDisTo(const IndividualType& indi, Problem *pro) {
			if (this->m_individuals.empty()) {
				return 0;
			}
			else {
				Real avg_dis(0);
				for (auto &it : this->m_individuals) {
					avg_dis += indi.variableDistance(*it, pro);
				}
				return avg_dis / static_cast<Real>(this->m_individuals.size());
			}
		}

		virtual void udpateCurRadius(Problem *pro) = 0;

		virtual void initRadius(Problem *pro) {
			udpateCurRadius(pro);
			m_initial_radius = m_cur_radius;
		}

		virtual bool judgeConverged() { return false; }

		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
			Real before_radius(m_cur_radius);
			int rf = TPopulation::evolve(pro, alg, rnd);
			udpateCurRadius(pro);
			if (m_cur_radius / (before_radius + 0.00001) < 0.9) m_count_radius_decre = 0;
			else m_count_radius_decre++;
			return rf;
		}
	};

}

#endif