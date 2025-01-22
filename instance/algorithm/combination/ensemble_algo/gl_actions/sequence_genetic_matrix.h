/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Framework of genetic learning (PopGL) algorithm
*
*********************************************************************************/

#ifndef SEQUENCE_GENETIC_MATRIX_H
#define SEQUENCE_GENETIC_MATRIX_H

#include"sequence_memory.h"
#include "../../../uncertianty/distance_calculator/fitness_weight_mapper.h"
#include "../../../../../utility/function/custom_function.h"


namespace ofec {

	class SequenceGeneticMatrixParameter :public SequenceMemoryParametersBase {
	public:
		enum class MatrixUpdateTag { kNone = 0, kMin };
		MatrixUpdateTag m_updateTag = MatrixUpdateTag::kMin;
	};

	template<typename TSequenceOperator>
	class SequenceGeneticMatrix :public SequenceMemory<TSequenceOperator> {
	protected:
		FitnessWeightMapper m_mapper;
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = SequenceGeneticMatrixParameter;

	public:
		//virtual void initialize(Problem *pro, Algorithm *alg, int id_param) {
		////	m_alpha = 0.9;
		//	SequenceMemory::initialize(pro, alg);
		//}

		virtual void initialize(const SequenceParametersBase& par) override {
			SequenceMemory::initialize(par);
			auto& cur_par = dynamic_cast<const ParType&>(par);
			m_matrix_update_tag = cur_par.m_updateTag;
		}


		virtual void globalUpdate(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) override{
			std::vector<ofec::Real> fitness(pop.size());
			for (int idx(0); idx < fitness.size(); ++idx) {
				fitness[idx] = pop[idx]->fitness();
			}
			m_mapper.reset();
			m_mapper.update(fitness);
			for (int idx(0); idx < fitness.size(); ++idx) {
				fitness[idx] = m_mapper.getPopFitness(fitness[idx]);
			}
			UTILITY::assignVVector<ofec::Real>(m_memory, 0);
			addWeight(pop, fitness);

			//for (int pidx(0); pidx < pop.size(); ++pidx) {
			//	auto& indi(*pop[pidx]);
			//	auto& weight(*vWeight[pidx]);
			//	for (auto& edge : indi.edges()) {
			//		m_memory[edge.first][edge.second] += weight;
			//	}
			//}

			if (m_matrix_update_tag == ParType::MatrixUpdateTag::kMin) {
				double minValue(fitness.front());
				for (auto& it : fitness) {
					minValue = std::min(it, minValue);
				}
				modifyMemoryMax(minValue);
			}
		}
	protected:
	//	double m_alpha = 0;
		ParType::MatrixUpdateTag m_matrix_update_tag = ParType::MatrixUpdateTag::kMin;
	};

}

#endif 