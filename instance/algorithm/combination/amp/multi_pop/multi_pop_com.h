/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/


#ifndef MULTI_POP_COM_H
#define MULTI_POP_COM_H

#include "gl_mp_pop.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"
#include "../../../template/combination/sequence/sequence_interpreter.h"
#include "../../../template/combination/sequence/sequence_operator.h"
#include "../../../../../core/algorithm/multi_population.h"
namespace ofec {
	template<typename TAdaptor>
	class MP_GL : public AlgorithmSeq {
	public:
		using AdaptorType = typename TAdaptor;
		using OpType = typename TAdaptor::OpType;
		using SolutionType = typename TAdaptor::SolutionType;
		using InterpreterType = typename TAdaptor::InterpreterType;
		using PopulationType = typename PopGLMpCom<TAdaptor>;
	protected:

		MultiPopulation<PopGLMpCom<TAdaptor>> m_pops;
		PopGLMpCom<TAdaptor> m_mainPop;
		Real m_InnerRadius = 0.8;

		WeightCalculator m_weight_calculator;
		DistanceCalculator<InterpreterType> m_distace_calculator;
		size_t m_iteration = 0;
		size_t m_maxEvolePop = 2;
		//	std::vector<std::unique_ptr<SolutionType>> m_history;
		//	std::vector<Solution> 
	public:
		virtual void initialize_()override;
		virtual void run_()override;
		// Í¨¹ý Algorithm ¼Ì³Ð
		virtual void record() override {}

		virtual void printPop(PopGLMpCom<TAdaptor>& pop);
	};

	template<typename TAdaptor>
	inline void MP_GL<TAdaptor>::initialize_()
	{
		AlgorithmSeq::initialize_();
		m_keep_candidates_updated = true;
		AlgorithmSeq::assignInterpreter<InterpreterType>();
		//		assignInterpreter<InterpreterType>(m_problem.get());
		assignOps<OpType>();
		m_distace_calculator.resize(getInterpreter<InterpreterType>());
		m_distace_calculator.initialize(this, m_problem.get());
		m_mainPop.assign(std::get<int>(m_param->at("population size")), m_problem.get());
		m_mainPop.resize(std::geimage sequence predictt<int>(m_param->at("population size")), m_problem.get());
		m_mainPop.initialize(m_problem.get(), this, m_random.get());
		m_mainPop.setRadiusDecrease(1.0);
		m_mainPop.setMaxIndiSize(1000);

	}

	template<typename TAdaptor>
	inline void MP_GL<TAdaptor>::run_()
	{
		auto& op = getOp<OpType>();
		auto& interpreter = getInterpreter<InterpreterType>();
		while (!terminating()) {
			for (auto& it : m_pops) it->evolve(m_problem.get(), this, m_random.get());
			for (int id(0); id < m_pops.size(); ++id) {
				m_pops[id].setId(id);				curInds[idx]->setId(-1);
			}
			//			curInds=std::move(m_mainPop.getIndis());
						//for (auto& it : curInds) {
						//	it.setActive(true);
						//	it.setId(-1);
						//}



			std::vector<std::unique_ptr<SolutionType>> newIndis;
			for (auto& popIter : m_pops) {
				newIndis.clear();
				for (auto& IndiIter : curInds) {
					if (IndiIter->isActive()) {
						auto dis2Pop =
							popIter->getCenter().variableDistance(*IndiIter, m_problem.get());
						if (dis2Pop <= m_InnerRadius) {
							IndiIter->setActive(false);
							if (dis2Pop <= popIter->getRadius()) {
								newIndis.emplace_back(std::move(IndiIter));
							}
						}
					}
				}
				popIter->mergePop(m_problem.get(), this, newIndis);
			}

			std::vector<int> popIdxs(m_pops.size());
			for (int popId(0); popId < popIdxs.size(); ++popId) {
				popIdxs[popId] = popId;
			}



			sort(popIdxs.begin(), popIdxs.end(), [&](int a, int b) {
				return op.better(m_problem.get(), this->at(a).getCenter(), this->at(b).getCenter());
			});

			for (auto& it : m_pops) {
				it->setActive(true);
			}

			for (int popId(0); popId < popIdxs.size(); ++popId)
			{
				auto& popA(m_pops[popIdxs[popId]]);
				if (popA != nullptr && popA->isActive()) {
					for (int popIj(popId + 1); popIj < popIdxs.size(); ++popIj) {
						auto& popB(m_pops[popIdxs[popIj]]);
						auto dis2Pop =
							popA->getCenter().variableDistance(popB->getCenter(), m_problem.get());
						if (dis2Pop < popA->getRadius()) {
							popA->mergePop(m_problem.get(), this, *popB);
							popB->setActive(false);
						}
						else if (dis2Pop < popB->getRadius()) {
							popB->setActive(false);
						}
					}
				}
			}
			//{
			//	int popId(0), curIdx(0);
			//	for ( curIdx=0; curIdx < m_pops.size(); ++curIdx) {
			//		while (popId < m_pops.size() && !m_pops[popId]->isActive()) {
			//			++popId;
			//		}
			//		if (popId == m_pops.size())break;
			//		if (popId != curIdx) {
			//			swap(m_pops[popId], m_pops[curIdx]);
			//		}
			//	}
			//	m_pops.resize(curIdx);
			//}

			{
				decltype(m_pops) newPop;
				for (auto& popIter : m_pops) {
					if (popIter->isActive()) {
						newPop.emplace_back(std::move(popIter));
					}
				}
				swap(m_pops, newPop);
			}

			{
				newIndis.clear();
				for (auto& it : curInds) {
					if (it->isActive()) {
						newIndis.emplace_back(std::move(it));
					}
				}
				swap(curInds, newIndis);
			}


			{
				sort(curInds.begin(), curInds.end(), [&](
					const std::unique_ptr<SolutionType>& a,
					const std::unique_ptr<SolutionType>& b) {
					return op.better(m_problem.get(), *a, *b);
				});

				int solId(0);
				for (; solId < curInds.size(); ++solId) {
					auto& solA(curInds[solId]);
					if (solA != nullptr) {
						newIndis.clear();
						for (int solIj(solId + 1); solIj < curInds.size(); ++solIj) {
							auto& solB(curInds[solIj]);
							if (solB != nullptr) {
								auto dis2Pop =
									solA->variableDistance(*solB, m_problem.get());
								if (dis2Pop < m_InnerRadius) {
									newIndis.emplace_back(std::move(solB));
								}
							}

						}
						newIndis.emplace_back(std::move(solA));
						m_pops.emplace_back(new PopulationType());
						m_pops.back()->initialize(m_problem.get(), this, m_random.get());
						m_pops.back()->assign(0, m_problem.get());
						m_pops.back()->setRadius(m_InnerRadius);
						m_pops.back()->mergePop(m_problem.get(), this, newIndis);
						m_pops.back()->addIndisAround(m_problem.get(), this);

						if (m_pops.size() > m_maxEvolePop) break;
					}

				}

			}


			{
				newIndis.clear();
				for (auto& popIter : m_pops) {
					newIndis.emplace_back(new SolutionType(popIter->getCenter()));
				}
				m_mainPop.resize(0, m_problem.get());
				m_mainPop.mergePop(m_problem.get(), this, newIndis);
				//	m_mainPop.resize(m_mainPop.size(), m_problem.get());

				++m_iteration;
				std::cout << m_iteration << std::endl;
				std::cout << "popSize\t" << m_pops.size() << std::endl;
				for (int id(0); id < m_pops.size(); ++id) {
					std::cout << "it\t" << id
						<< "\tbestOjb\t" << m_pops[id]->getCenter().objective()[0]
						<< "\tdis\t" << m_pops[id]->getRadius() << std::endl;
				}
			}



#ifdef OFEC_DEMO
			updateBuffer();
#endif
			//	m_problem->showInfomations(this);
		}

	}



	template<typename TAdaptor>
	void MP_GL<TAdaptor>::printPop(PopGLMpCom<TAdaptor>& m_pop) {
		Real minObj(m_pop.begin()->get()->objective()[0]);
		Real maxObj(minObj);
		for (auto& it(m_pop.begin()); it != m_pop.end(); ++it) {
			minObj = std::min((*it)->objective()[0], minObj);
			maxObj = std::max((*it)->objective()[0], maxObj);
		}
		std::cout << "minOjb\t" << minObj << "\tmaxOjb\t" << maxObj << std::endl;

		{
			std::vector<Real> weight;
			{
				std::vector<SolutionBase*> pop(m_pop.size());
				for (int idx(0); idx < pop.size(); ++idx) {
					pop[idx] = &m_pop[idx];
				}
				m_weight_calculator.calWeight(m_problem.get(), pop, weight);
			}

			{
				std::vector<SolutionType*> pop(m_pop.size());
				for (int idx(0); idx < pop.size(); ++idx) {
					pop[idx] = &m_pop[idx];
				}
				m_distace_calculator.addDistance(m_problem.get(), getInterpreter<InterpreterType>(), pop, weight);
				m_distace_calculator.printInfo();
			}
			std::cout << "numImprove\t" << m_pop.numImprove() << std::endl;
		}
		int pidx(0);
		for (int idx(0); idx < m_pop.size(); ++idx) {
			if (m_pop[idx].dominate(m_pop[pidx], m_problem.get())) {
				pidx = idx;
			}
		}
		std::cout << "solution id\t" << pidx << std::endl;
		m_problem->printfSolution(this, m_pop[pidx]);
	}

}

#endif 