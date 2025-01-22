/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/


#ifndef GL_SEQ_H
#define GL_SEQ_H


#include"gl_pop_seq.h"
#include "../../../../core/problem/uncertainty/dynamic.h"

#include"../../template/combination/sequence/sequence_algorithm.h"
#include"../../template/combination/sequence/sequence_interpreter.h"
#include"../../template/combination/sequence/sequence_operator.h"
#include "../../../../instance/problem/combination/selection_problem/selection_problem.h"
#include"../../template/combination/multi_population/distance_calculator.h"
#include "../../combination/ensemble_algo/gl_actions/weight_calculator.h"

#include "../ensemble_algo/gl/adaptor_gl_seq.h"

namespace ofec {
	template<typename TAdaptor>
	class GL_Seq : public AlgorithmSeq {
	public:
		using AdaptorType = typename TAdaptor;
		using OpType = typename TAdaptor::OpType;
		using SolutionType = typename TAdaptor::SolutionType;
		using InterpreterType = typename TAdaptor::InterpreterType;
	protected:
		PopGLSeq<AdaptorType> m_pop;
		WeightCalculator m_weight_calculator;
		DistanceCalculator<InterpreterType> m_distace_calculator;
		Real m_radius = 1.0;
	public:
		virtual void initialize_()override;
		virtual void run_()override;
		// Í¨¹ý Algorithm ¼Ì³Ð
		virtual void record() override {
			GET_DYN_SP(m_problem.get())->recordSol(*candidates().front());
		//	virtual void record() = 0;
		}



	};

	template<typename TAdaptor>
	inline void GL_Seq<TAdaptor>::initialize_()
	{
		AlgorithmSeq::initialize_();
		m_keep_candidates_updated = true;
		assignInterpreter<InterpreterType>();
		assignOps<OpType>();
		m_distace_calculator.resize(getInterpreter<InterpreterType>());
		m_distace_calculator.initialize(m_problem.get(), this);
		m_pop.assign(std::get<int>(m_param->at("population size")),m_problem.get());
		m_pop.resize(std::get<int>(m_param->at("population size")), m_problem.get());
		m_pop.initialize(m_problem.get(), this, m_random.get());
		

	}

	template<typename TAdaptor>
	inline void GL_Seq<TAdaptor>::run_()
	{
		while (!terminating()) {
		//	std::cout << "alg id\t" << this << std::endl;
			if (m_pop.getRadius() < 0.1) {
				m_pop.setRadius(1.0);
			}
			m_pop.evolve(m_problem.get(), this, m_random.get());
			//m_pop.updateBest(m_problem.get());
		//	std::cout << evaluations() << std::endl;
			m_problem->showInfomations(this);


		//
//			m_problem->printfSolution(this, m_pop.best().front());
			
			//Real minObj(m_pop.begin()->get()->objective()[0]);
			//Real maxObj(minObj);
			//for (auto& it(m_pop.begin()); it != m_pop.end(); ++it) {
			//	minObj = std::min((*it)->objective()[0], minObj);
			//	maxObj = std::max((*it)->objective()[0], maxObj);
			//}
			//std::cout << "minOjb\t" << minObj << "\tmaxOjb\t" << maxObj << std::endl;
	/*		if (m_problem->hasTag(ProblemTag::kDOP)) {
				GET_DOP(m_problem.get())->setCurSolution(m_pop.getCenter());
			}*/
			//m_problem->setCurSolution(m_pop.getCenter());
			
			//{
			//	std::vector<Real> weight;
			//	{
			//		std::vector<SolutionBase*> pop(m_pop.size());
			//		for (int idx(0); idx < pop.size(); ++idx) {
			//			pop[idx] = &m_pop[idx];
			//		}
			//		m_weight_calculator.calWeight(m_problem.get(), pop, weight);
			//	}

			//	{
			//		std::vector<SolutionType*> pop(m_pop.size());
			//		for (int idx(0); idx < pop.size(); ++idx) {
			//			pop[idx] = &m_pop[idx];
			//		}
			//		m_distace_calculator.addDistance(m_problem.get(), getInterpreter<InterpreterType>(), pop, weight);
			//		m_distace_calculator.printInfo();
			//	}
			//	std::cout <<"numImprove\t" << m_pop.numImprove() << std::endl;
			//}


			//int pidx(0);
			//for (int idx(0); idx < m_pop.size(); ++idx) {
			//	if (m_pop[idx].dominate(m_pop[pidx],m_problem.get())){
			//		pidx = idx;
			//	}
			//}
			//std::cout << "solution id\t" << pidx << std::endl;
			//m_problem->printfSolution(this, m_pop[pidx]);
			//m_pop.setRadius(m_radius);
			//m_pop.shrinkPop(m_problem.get(), this, m_random.get());
			//m_radius -= 0.01;
			//if (m_radius <0.1) {
			//	m_radius = 0.1;
			//}

		//	m_problem->showInfomations(this);
		}
	}

}

#endif 