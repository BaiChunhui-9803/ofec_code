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

#ifndef SEQUENCE_FREQUENCY_H
#define SEQUENCE_FREQUENCY_H

#include"sequence_memory.h"

namespace ofec {
	template<typename TSequenceOperator>
	class SequenceFrequency :public SequenceMemory<TSequenceOperator> {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
	    using ParType = typename  SequenceMemory<TSequenceOperator>::ParType;
	public:
		virtual void initialize(const SequenceParametersBase& par) override {
			SequenceMemory::initialize(par);
			//m_initVal = cur_par.m_initVal;
			//for (auto& it : m_memory) {
			//	std::fill(it.begin(), it.end(), m_initVal);
			//}
		}

		virtual void globalUpdate(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) override {
			for (auto& it : m_memory) {
				for (auto& it2 : it) {
					it2 = it2 * m_alpha + (1.0 - m_alpha);
				}
			}
			double weight = (1.0 / pop.size()) * (1.0 - m_alpha);
			for (int pidx(0); pidx < pop.size(); ++pidx) {
				auto& indi(*pop[pidx]);
				for (auto& edge : indi.edges()) {
					m_memory[edge.first][edge.second] -= weight;
				}
			}
		
			//updateMaxPro();
		}
	protected:
		double m_alpha = 0;
	};

}

#endif