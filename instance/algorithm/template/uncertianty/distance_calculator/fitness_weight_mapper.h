
#ifndef FITNESS_WEIGHT_MAPPER_H
#define FITNESS_WEIGHT_MAPPER_H

#include<vector>
#include<list>
#include<iomanip>

namespace ofec {
	class FitnessWeightMapper {
	public:
		enum class mode { singleObjSigmoid, singleNormalize, origin };
		mode m_mode = mode::singleObjSigmoid;
	protected:
		double m_max_fitness = 0, m_min_fitness = 0, m_gap = 0;

		double m_max_weight = 0;
		double m_min_weight = 0;
	public:

		double getMaxWeight()const {
			return m_max_weight;
		}
		double getMinWeight()const {
			return m_min_weight;
		}
		virtual void reset() {
			m_max_fitness = -std::numeric_limits<double>::max();
			m_min_fitness = -m_max_fitness;
			m_gap = std::numeric_limits<double>::max();
		}

		virtual void update(const std::vector<double>& fitness) {
			for (auto& fit : fitness) {
				if (m_max_fitness < fit)
					m_max_fitness = fit;
				if (m_min_fitness > fit)
					m_min_fitness = fit;
			}
			m_gap = m_max_fitness - m_min_fitness + 1e-5;

			m_max_weight = getPopFitness(m_max_fitness);
			m_min_weight = getPopFitness(m_min_fitness);
			
		}

		template<typename TIndi>
		void update(std::vector<std::unique_ptr<TIndi>>& pop) {
			for (int i = 1; i < pop.size(); ++i) {
				if (m_max_fitness < pop[i]->fitness())
					m_max_fitness = pop[i]->fitness();
				if (m_min_fitness > pop[i]->fitness())
					m_min_fitness = pop[i]->fitness();
			}
			Real gap = m_max_fitness - m_min_fitness + 1e-5;
		}

		double getPopFitness(double fitness) const{
			if (fitness > m_max_fitness) fitness = m_max_fitness;
			else if (fitness < m_min_fitness) fitness = m_min_fitness;
			double fitVal(0);
			switch (m_mode) {
			case mode::singleObjSigmoid: {
				fitVal = (fitness - m_min_fitness + 1e-5) / m_gap;
				fitVal = 1 / (1 + exp(-fitVal));
			}
									   break;
			case mode::singleNormalize: {
				fitVal = (fitness - m_min_fitness + 1e-5) / m_gap;
			}
									  break;
			case mode::origin: {
				fitVal = fitness;
			}
			default: break;
			}
			return fitVal;
		}

	};
}


#endif