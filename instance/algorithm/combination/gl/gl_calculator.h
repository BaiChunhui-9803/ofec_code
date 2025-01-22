#ifndef GL_SEQ_CALCULATOR_H
#define GL_SEQ_CALCULATOR_H

#include"../../../../../core/definition.h"
#include "../../../../../core/problem/encoding.h"
#include<vector>
#include<list>

namespace ofec {
	class GLSeqCalculator {
	public:
		enum class WeightMode { singleObjSigmoid, singleNormalize };

		void calWeight(Problem *pro, const std::vector<SolutionBase*>& pop,
			std::vector<Real>& weight, 
			WeightMode mode= WeightMode::singleObjSigmoid,
			int objIdx=0) {
			weight.resize(pop.size());
			if (pop.empty())return;
			switch (mode)
			{
			case weightMode::singleObjSigmoid:
			{
				Real m_preMemoryMaxObj = pop.front()->objective(objIdx);
				Real m_preMemoryMinObj = pop.front()->objective(objIdx);

				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->objective(objIdx))
						m_preMemoryMaxObj = pop[i]->objective(objIdx);
					if (m_preMemoryMinObj > pop[i]->objective(objIdx))
						m_preMemoryMinObj = pop[i]->objective(objIdx);
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
					if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
						weight[i] = (m_preMemoryMaxObj - pop[i]->objective(objIdx) + 1e-5) / gap;
					else
						weight[i] = (pop[i]->objective(objIdx) - m_preMemoryMinObj + 1e-5) / gap;
					weight[i] = 1 / (1 + exp(-weight[i]));
				}

			}
			break;
			case WeightMode::singleNormalize:
			{
				Real m_preMemoryMaxObj = pop.front()->objective(objIdx);
				Real m_preMemoryMinObj = pop.front()->objective(objIdx);

				for (int i = 1; i < pop.size(); ++i) {
					if (m_preMemoryMaxObj < pop[i]->objective(objIdx))
						m_preMemoryMaxObj = pop[i]->objective(objIdx);
					if (m_preMemoryMinObj > pop[i]->objective(objIdx))
						m_preMemoryMinObj = pop[i]->objective(objIdx);
				}
				Real gap = m_preMemoryMaxObj - m_preMemoryMinObj + 1e-5;
				for (int i = 0; i < pop.size(); ++i) {
					if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
						weight[i] = (m_preMemoryMaxObj - pop[i]->objective(objIdx) + 1e-5) / gap;
					else
						weight[i] = (pop[i]->objective(objIdx) - m_preMemoryMinObj + 1e-5) / gap;
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