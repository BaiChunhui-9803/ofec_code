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

#ifndef SEQUENCE_MEMORY_H
#define SEQUENCE_MEMORY_H
#include"../../../../../core/definition.h"
#include "../../../template/combination/sequence/sequence_action.h"

#include<vector>
#include<memory>

namespace ofec {

	class SequenceMemoryParametersBase :public SequenceParametersBase {
	public:
		int m_initVal = 1.0;
		int m_alpha = 0.9;
	};


	template<typename TSequenceOperator>
	class SequenceMemory :public SequenceAction <TSequenceOperator>{
	public:
		using OpType = typename TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = typename SequenceMemoryParametersBase;
		
	protected:
		int this = -1;
		int m_problem.get() = -1;

	protected:

		void modifyMemoryTimes(Real rho = 0.0) {

			for (auto& it : m_memory) {
				for (auto& it2 : it) {
					it2 *= rho;
				}
			}
		}
		void modifyMemoryMax(Real val = 0.0) {

			for (auto& it : m_memory) {
				for (auto& it2 : it) {
					it2 = std::max(it2, val);
				}
			}
		}

		void modifyMemoryMin(Real val = 0.0) {

			for (auto& it : m_memory) {
				for (auto& it2 : it) {
					it2 = std::min(it2, val);
				}
			}
		}

		inline void addWeight(
			const std::vector<SolutionType *>& pop, 
			const std::vector<Real>& weight)
		{

			for (int pidx(0); pidx < weight.size(); ++pidx) {
				auto& indi(*pop[pidx]);
				for (auto& edge : indi.edges()) {
					m_memory[edge.first][edge.second] += weight[pidx];
				}
			}
		}

		//virtual void updateMaxPro() {
		//	m_max_pro.resize(mat_size.size());
		//	double tot_pro(0), max_pro(0);
		//	for (int idx(0); idx < m_max_pro.size(); ++idx) {
		//		auto& it(m_memory[idx]);
		//		tot_pro = 0;
		//		max_pro = 0;
		//		for (auto& it2 : it) {
		//			tot_pro += it2;
		//			max_pro = std::max(max_pro, it2);
		//		}
		//		if (tot_pro > 0) {
		//			m_max_pro[idx] = max_pro / tot_pro;
		//		}
		//		else {
		//			m_max_pro[idx] = 0;
		//		}
		//	}
		//}
	public:
		SequenceMemory() = default;
		~SequenceMemory() = default;

		virtual void initialize(const SequenceParametersBase& par) override {
			
			m_flagStepLearn = true;
			auto& cur_par = dynamic_cast<const ParType&>(par);
			this = cur_par.this;
			m_problem.get() = cur_par.m_problem.get();
			m_initVal = cur_par.m_initVal;
			auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
			resize(interpreter.getMatrixSize());
		}

		virtual void resize(const std::vector<int>& mat_size, double* initVal = nullptr)override {
			if (initVal != nullptr) {
				m_initVal = *initVal;
			}
			m_memory.resize(mat_size.size());
			for (int idx(0); idx < mat_size.size(); ++idx) {
				m_memory[idx].resize(mat_size[idx], m_initVal);
			}
		//	updateMaxPro();
		}


		virtual int learn(
			SolutionType& cur,
			 Problem *pro, Random *rnd,
			const InterpreterType& interpreter,
			const OpType& op,
			const std::vector<SolutionType*>& pop,
			int id_indiv
		)const override {
			std::function<Real(const SolutionType& cur, int to)>
				pro_fun = [&](const SolutionType& cur, int to) {
				return getPro(pro, interpreter, cur, to);
			};
			//= std::bind(&SequenceMemory::getPro, this,
			//	pro, interpreter, std::placeholders::_1, std::placeholders::_2);
			return op.learn_from_global(cur, rnd, pro, interpreter, pro_fun);

		}


		virtual Real getPro(Problem *pro,
			const InterpreterType& interpreter,
			const SolutionType& cur, int to) const{
			return m_memory[interpreter.curPositionIdx(pro, cur)][to];
		}
		//virtual Real getMaxPro(Problem *pro,
		//	const InterpreterType& interpreter,
		//	const SolutionType& cur)override {
		//	return m_max_pro[interpreter.curPositionIdx(pro, cur)];
		//}
	protected:
		double m_initVal = 1.0;
		std::vector<std::vector<ofec::Real>> m_memory;
		
	//	std::vector<double> m_max_pro;
		//std::function<Real(const SolutionType& cur, int to)> m_pro_fun;
	};

}

#endif