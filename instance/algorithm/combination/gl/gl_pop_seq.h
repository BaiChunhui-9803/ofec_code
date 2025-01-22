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

#ifndef OFEC_GL_SEQ_H
#define OFEC_GL_SEQ_H

#include  "../../template/framework/gl/gl_pop.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	template<typename TAdaptor>
	class PopGLSeq :virtual public PopGL<typename TAdaptor::SolutionType> {
	public:
		using AdaptorType = typename TAdaptor;
		using OpType = typename TAdaptor::OpType;
		using SolutionType = typename TAdaptor::SolutionType;
		using InterpreterType = typename TAdaptor::InterpreterType;

	public:
        void initialize(Problem *pro, Algorithm *alg, Random *rnd);
		virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		virtual void setRadius(Real radius) {
			dynamic_cast<AdaptorType&>(*m_adaptor).setRadius(radius);
		}

		AdaptorType& getAdaptor() {
			return dynamic_cast<AdaptorType&>(*m_adaptor);
		}
		const SolutionType& getCenter() {
			return dynamic_cast<AdaptorType&>(*m_adaptor).getCenter();
		}
		
		virtual Real getRadius() {
			return dynamic_cast<AdaptorType&>(*m_adaptor).getRadius();
		}
		virtual void shrinkPop(Problem *pro, Algorithm *alg, Random *rnd);

		virtual Real popRadius(Problem *pro) {
			auto& center(dynamic_cast<AdaptorType&>(*m_adaptor).getCenter());
			 Real curPosAvgDis = 0;
			for (auto& popIter : this->m_individuals) {
				curPosAvgDis += center.variableDistance(*popIter, pro);
			}
			curPosAvgDis /= m_individuals.size();
			return curPosAvgDis;
		}
		void setRadiusDecrease(Real val) {
			m_radius_decrease = val;
		}
	protected:
	//	int m_numImp = 0;
		Real m_radius_decrease = 0.9;
		int m_NotConvergeTime = 0;
		int m_NotConvergeTimeThread = 50;
		bool m_flag_shrink_circle = false;

	};
	template<typename TSequenceOperator>
	void PopGLSeq<TSequenceOperator>::shrinkPop(Problem *pro, Algorithm *alg, Random *rnd){

		dynamic_cast<AdaptorType&>(*m_adaptor).shrinkPop(pro, alg, rnd, this->m_individuals);
		for (auto& i : m_individuals)
			updateBest(*i, pro);
	}


	template<typename TSequenceOperator>
	void PopGLSeq<TSequenceOperator>::initialize(Problem *pro, Algorithm *alg , Random *rnd) {

		m_ms = UpdateScheme::ci;
		m_adaptor.reset(new AdaptorType());
		Population<SolutionType>::initialize(pro, rnd);
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<OpType>());
		for (auto& it : m_individuals) {
			it->reset();
			interpreter.stepFinal(pro, *it);
			op.evaluate(pro, alg, rnd, *it);
		}
		getAdaptor().initialize(pro, alg);
		dynamic_cast<AdaptorType&>(*m_adaptor).updateCenter(pro, alg, m_individuals);
		initializeMemory(pro,alg);
		for (int i = 0; i < this->size(); i++) {
			m_offspring.push_back(*this->m_individuals[i]);
		}
	}
	template<typename TSequenceOperator>
	int PopGLSeq<TSequenceOperator>::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		updateMemory(pro, alg, true);
		m_adaptor->createSolution(pro,alg, rnd,m_individuals, m_offspring);
		int rf = update(pro, alg);
		for (auto& i : m_individuals)
				updateBest(*i, pro);

		m_iteration++;

		updateMemory(pro, alg, true);
		//getAdaptor().updatePointDensity(pro, alg, m_individuals);

#ifdef OFEC_DEMO

		auto& inter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		for (auto& it : m_individuals) inter.updateEdges(pro, *it);
		getAdaptor().updateDensityMatrix(m_individuals);
		getAdaptor().updatePointDensity(pro, alg, m_individuals);
		getAdaptor().updateOptPointInfo(pro, alg);
		//ofec_demo::g_buffer->appendAlgBuffer(alg);
#endif


		if (m_flag_shrink_circle) {
			Real curRadius = popRadius(pro);
			Real setRadiusValue = getRadius() * m_radius_decrease;

			//	std::cout << "curRadius\t" << curRadius << "\tsetRadius\t" << setRadiusValue << std::endl;
			if (curRadius <= setRadiusValue) {
				m_NotConvergeTime = 0;
				setRadius(setRadiusValue);
			}
			else {
				++m_NotConvergeTime;
				if (m_NotConvergeTime >= m_NotConvergeTimeThread) {
					m_NotConvergeTime = 0;
					curRadius = setRadiusValue;
					setRadius(curRadius);
					shrinkPop(pro, alg, rnd);
					if (m_NotConvergeTime == 0) {
						//	std::cout << "shrinkPop\t" << std::endl;

					}
				}
			}

		}
		return rf;
	}


	/*

	*/

}






#endif

