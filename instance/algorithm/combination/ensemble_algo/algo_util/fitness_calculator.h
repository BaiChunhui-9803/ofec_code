
#ifndef Fitness_CALCULATOR_H
#define Fitness_CALCULATOR_H

#include"../../../../../core/definition.h"
#include "../../../../../core/problem/encoding.h"

#include<vector>
#include<list>

namespace ofec {
	class FitnessCalculator {
	public:
		enum class mode { singleObjSigmoid, singleNormalize, origin};
		mode m_mode = mode::singleObjSigmoid;

	protected:
		double m_max_fitness = 0, m_min_fitness = 0,m_gap = 0;
		int m_problem.get() = -1;
	public:

		virtual void initialize(Problem *pro) {
			m_problem.get() = pro;
			reset();
		}

		virtual void reset() {
			m_max_fitness = -std::numeric_limits<double>::max();
			m_min_fitness = -m_max_fitness;
			m_gap = std::numeric_limits<double>::max();
		}

		virtual double ObjectiveToFitness(const SolutionBase& sol) {
			if (m_problem->optimizeMode(0) == OptimizeMode::kMinimize) {
				return -sol.objective()[0];
			}
			else return sol.objective()[0];
		}


		virtual void update(const std::vector<double>& fitness) {
			for (auto& fit : fitness) {
				if (m_max_fitness < fit)
					m_max_fitness = fit;
				if (m_min_fitness > fit)
					m_min_fitness = fit;
			}
			m_gap = m_max_fitness - m_min_fitness + 1e-5;
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

		double getPopFitness(double fitness) {
			if (fitness > m_max_fitness) fitness = m_max_fitness;
			else if (fitness < m_min_fitness) fitness = m_min_fitness;
			double fitVal (0);
			switch (m_mode){
			case mode::singleObjSigmoid:{
				fitVal = (fitness - m_min_fitness + 1e-5) / m_gap;
				fitVal = 1 / (1 + exp(-fitVal));
			}
			break;
			case mode::singleNormalize:{
				fitVal = (fitness - m_min_fitness + 1e-5) / m_gap;
			}
			break;
			case mode::origin: {
				fitVal = fitness;
			}
			default: break;}
			return fitVal;
		}

	};
}


#endif