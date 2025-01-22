
#ifndef WEIGHT_CALCULATOR_H
#define WEIGHT_CALCULATOR_H

#include"../../../../../core/definition.h"
#include "../../../../../core/problem/encoding.h"

#include<vector>
#include<list>

namespace ofec {
	class WeightCalculator {
		enum class mode { singleObjSigmoid , singleNormalize};
		mode m_mode = mode::singleObjSigmoid;
	public:

		void calWeight(Problem *pro,const std::vector<SolutionBase*>& pop,
			std::vector<Real>& weight) {
			weight.resize(pop.size());
			if (pop.empty())return;
			switch (m_mode)
			{
			case WeightCalculator::mode::singleObjSigmoid:
			{
				Real m_preMemoryMaxObj = pop.front()->objective(0);
				Real m_preMemoryMinObj = pop.front()->objective(0);
				
				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->objective(0))
						m_preMemoryMaxObj = pop[i]->objective(0);
					if (m_preMemoryMinObj > pop[i]->objective(0))
						m_preMemoryMinObj = pop[i]->objective(0);
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
					if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
						weight[i] = (m_preMemoryMaxObj - pop[i]->objective(0) + 1e-5) / gap;
					else
						weight[i] = (pop[i]->objective(0) - m_preMemoryMinObj + 1e-5) / gap;
					weight[i] = 1 / (1 + exp(-weight[i]));
				}

			}
				break;
			case WeightCalculator::mode::singleNormalize:
			{
				Real m_preMemoryMaxObj = pop.front()->objective(0);
				Real m_preMemoryMinObj = pop.front()->objective(0);

				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->objective(0))
						m_preMemoryMaxObj = pop[i]->objective(0);
					if (m_preMemoryMinObj > pop[i]->objective(0))
						m_preMemoryMinObj = pop[i]->objective(0);
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
					if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
						weight[i] = (m_preMemoryMaxObj - pop[i]->objective(0) + 1e-5) / gap;
					else
						weight[i] = (pop[i]->objective(0) - m_preMemoryMinObj + 1e-5) / gap;
				}

			}
				break;
			default:
				break;
			}


		}
	
	
		void calWeightFitness(Problem *pro, const std::vector<SolutionBase*>& pop,
			std::vector<Real>& weight) {
			weight.resize(pop.size());
			if (pop.empty())return;
			switch (m_mode)
			{
			case WeightCalculator::mode::singleObjSigmoid:
			{
				Real m_preMemoryMaxObj = pop.front()->fitness();
				Real m_preMemoryMinObj = pop.front()->fitness();

				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->fitness())
						m_preMemoryMaxObj = pop[i]->fitness();
					if (m_preMemoryMinObj > pop[i]->fitness())
						m_preMemoryMinObj = pop[i]->fitness();
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
					weight[i] = (pop[i]->fitness() -m_preMemoryMinObj + 1e-5) / gap;
					weight[i] = 1 / (1 + exp(-weight[i]));
				}

			}
			break;
			case WeightCalculator::mode::singleNormalize:
			{
				Real m_preMemoryMaxObj = pop.front()->fitness();
				Real m_preMemoryMinObj = pop.front()->fitness();

				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->fitness())
						m_preMemoryMaxObj = pop[i]->fitness();
					if (m_preMemoryMinObj > pop[i]->fitness())
						m_preMemoryMinObj = pop[i]->fitness();
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
						weight[i] = (pop[i]->fitness() - m_preMemoryMinObj + 1e-5) / gap;
				}

			}
			break;
			default:
				break;
			}


		}
	
	};
}


#endif