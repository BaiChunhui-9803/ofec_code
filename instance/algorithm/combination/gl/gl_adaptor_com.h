/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Changhe Li & Yiya Diao
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/ofec for more information
*
*-------------------------------------------------------------------------------
* Operator adaptator of genetic learning algorithm
*
*********************************************************************************/

#ifndef OFEC_GL_ADAPTER_SEQ_H
#define OFEC_GL_ADAPTER_SEQ_H


#include"../../template/framework/gl/gl_adaptor.h"
#include"../../template/combination/sequence/sequence_algorithm.h"
#include "../../../../utility/functional.h"




namespace ofec {
	template<typename TSequenceOperator>
	class AdaptorGLSeq : public AdaptorGL<typename TSequenceOperator::SolutionType> {

	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;

	protected:
		std::vector<std::vector<Real>> m_pro;
		enum class proUpdate { origin, min };
		proUpdate m_pro_update = proUpdate::min;
		int m_localSearchTimes = 30;
		Real m_radius = 1.0;
		SolutionType m_center;
		bool m_localSearch = true;
		bool m_dynamicHandle = true;


		// for test

	//	std::list<std::unique_ptr<SolutionType>> m_discard_indis;
	//	std::list<Real> opVar;
	protected:
		void resetPro(const InterpreterType& interpreter) {
			auto& mat(interpreter.getMatrixSize());
			m_pro.resize(mat.size());
			for (int idx(0); idx < mat.size(); ++idx) {
				m_pro[idx].resize(mat[idx]);
				std::fill(m_pro[idx].begin(), m_pro[idx].end(), 0);
			}
		}

		void modifyWeightMin(const std::vector<Real>& weight) {
			Real minValue(weight.front());
			for (auto& it : weight) minValue = std::min(minValue, it);
			for (auto& it : m_pro) {
				for (auto& it2 : it) {
					it2 = std::max(it2, minValue);
				}
			}
		}
		//std::list<std::unique_ptr<SolutionType>>& getDiscardIndis() {
		//	return m_discard_indis;
		//}
	public:
		AdaptorGLSeq() {
			m_alpha = 0.8;
		}
		virtual ~AdaptorGLSeq() {}


		//double getPro(SolutionType& indi) {
		//	indi.calEdges();
		//}

		virtual void initialze() {
			//		m_discard_indis.clear();
		}
		const SolutionType& getCenter()const {
			return m_center;
		}


		SolutionType& getCenter() {
			return m_center;
		}
		void setCenter(const SolutionType& center) {
			m_center = center;
		}
		Real getRadius() {
			return m_radius;
		}
		virtual void setRadius(Real radius) {
			m_radius = radius;
		}
		virtual void updateProbability(Problem *pro, Algorithm *alg,
			std::vector<std::unique_ptr<SolutionType>>& pop,
			const std::vector<Real>& weight,
			const std::vector<int>* index = nullptr)override;
		virtual void updateProbability(Problem *pro, Algorithm *alg,
			std::vector<SolutionType>& pop,
			const std::vector<Real>& weight,
			const std::vector<int>* index = nullptr)override;
		virtual void createSolution(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<std::unique_ptr<SolutionType>>& pop,
			std::vector<SolutionType>& offspring)override;
		virtual int updateSolution(Problem *pro, Algorithm *alg,
			std::vector<std::unique_ptr<SolutionType>>& pop,
			std::vector<SolutionType>& offspring, int& num_improve)override;


		virtual void shrinkPop(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop);
		virtual void updateCenter(Problem *pro, Algorithm *alg, std::vector<std::unique_ptr<SolutionType>>& pop) {
			auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
			int bestId = calIdx<std::unique_ptr<SolutionType>>(pop, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
				return op.better(pro, *a, *b);
			});
			if (bestId == -1)return;
			m_center = *pop[bestId];
			//if (!opVar.empty()) {
			//	if (opVar.back() - m_center.objective()[0] < -0.01) {
			//		std::cout << "center\t" << m_center.objective()[0] << std::endl;
			//		int stop = -1;
			//	}
			//}
			//opVar.push_back(m_center.objective()[0]);

		}
	};


	template<typename TSequenceOperator>
	inline void AdaptorGLSeq<TSequenceOperator>::shrinkPop(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop)
	{
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
		for (auto& it : pop) {
			if (op.learn_from_other(*it, rnd, pro, interpreter, m_center,m_radius)) {
			//	std::cout << "distance2center\t" << m_center.variableDistance(*it,pro) << std::endl;
				it->reset();
				interpreter.stepFinal(pro, *it);
				op.evaluate(pro, alg, rnd, *it);
			}
		}
		//updateCenter(pro, pop);
	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeq<TSequenceOperator>::updateProbability
	(Problem *pro, Algorithm *alg, 
		std::vector<std::unique_ptr<SolutionType>>& pop, 
		const std::vector<Real>& weight, 
		const std::vector<int>* index)
	{	
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		resetPro(interpreter);
		if (index != nullptr) {
			for (auto& pidx : *index) {
				auto& indi(pop[pidx]);
				interpreter.updateEdges(pro, *indi);
				for (auto& edge : indi->edges()) {
					m_pro[edge.first][edge.second] += weight[pidx];
				}
			}
		}
		else {
			for (int pidx(0); pidx < weight.size();++pidx) {
				auto& indi(pop[pidx]);
				interpreter.updateEdges(pro, *indi);
				for (auto& edge : indi->edges()) {
					m_pro[edge.first][edge.second] += weight[pidx];
				}
			}
			
		}
		if (m_pro_update == proUpdate::min) {
			modifyWeightMin(weight);
		}
		//auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());

	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeq<TSequenceOperator>::updateProbability(Problem *pro, Algorithm *alg, std::vector<SolutionType>& pop, const std::vector<Real>& weight, const std::vector<int>* index)
	{

		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		resetPro(interpreter);
		if (index != nullptr) {
			for (auto& pidx : *index) {
				auto& indi(pop[pidx]);
				interpreter.updateEdges(pro, indi);
				for (auto& edge : indi.edges()) {
					m_pro[edge.first][edge.second] += weight[pidx];
				}
			}
		}
		else {
			for (int pidx(0); pidx < weight.size(); ++pidx) {
				auto& indi(pop[pidx]);
				interpreter.updateEdges(pro, indi);
				for (auto& edge : indi.edges()) {
					m_pro[edge.first][edge.second] += weight[pidx];
				}
			}

		}
		if (m_pro_update == proUpdate::min) {
			modifyWeightMin(weight);
		}
		//auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());

	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeq<TSequenceOperator>::createSolution(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop, std::vector<SolutionType>& offspring)
	{
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());

		for(int  offId(0);offId< offspring.size();++offId)
		{
			auto& it = offspring[offId];
			it.reset();
			interpreter.stepInit(pro, it);
			while (!interpreter.stepFinish(pro, it)) {
				int next(-1);
				if (interpreter.stepFeasible(pro, it)) {
					if (rnd->uniform.next() < m_alpha
						/*std::max(m_alpha,1.0-m_radius/it->variable().size())*/) {
						next = op.learn_from_other(it, rnd, pro, interpreter,*pop[offId]);
					}
					if (next == -1) {
						auto weightFun = [&](const SolutionType& indi, int to) {
							return m_pro[interpreter.curPositionIdx(pro, indi)][to];
						};
						next = op.learn_from_global(it, rnd, pro, interpreter, weightFun);
					}
				}
				if (next != -1) {
					interpreter.stepNext(pro, it, next);
				}
				else {
					interpreter.stepBack(pro, it);
				}
			}
			op.learn_from_other(it, rnd, pro, interpreter, m_center, m_radius);
			interpreter.stepFinal(pro, it);
		}
		
	}
	template<typename TSequenceOperator>
	inline int AdaptorGLSeq<TSequenceOperator>::updateSolution(Problem *pro, Algorithm *alg, std::vector<std::unique_ptr<SolutionType>>& pop, std::vector<SolutionType>& offspring, int& num_improve)
	{

		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
		auto& DNStrategy = GET_ASeq(alg).getDNstrategy();
		int rf = 0;
		if (m_dynamicHandle) {
			std::vector<Real> origin_objs;
			for (int pidx(0); pidx < pop.size(); ++pidx) {
				origin_objs = pop[pidx]->objective();
				rf |= op.evaluate(pro, alg, alg.idRandom(), *pop[pidx]);
				{
					if (DNStrategy.judgeChange(pro, origin_objs, pop[pidx]->objective())) {
						pop[pidx]->setLocalSearch(false);
					}
				}
			//	std::cout << "obj different\t" << origin_objs.front() - pop[pidx]->objective().front() << std::endl;
			}
		}
		for (int pidx(0); pidx < pop.size(); ++pidx) {
		//	if (offspring[pidx].same(*pop[pidx], pro))continue;
			rf |= op.evaluate(pro, alg, alg.idRandom(), offspring[pidx]);
		}
		// getBestIndex
		{

			int popBestIndex = calIdx<std::unique_ptr<SolutionType>>(pop, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
				return op.better(pro, *a, *b);
			});
			int offBestIndex = calIdx<SolutionType>(offspring, 
				[&](const SolutionType& a, const SolutionType& b) {
				return op.better(pro, a, b);
			});

			if (op.better(pro, *pop[popBestIndex], offspring[offBestIndex])) {
				m_center = *pop[popBestIndex];
			}
			else {
				m_center = offspring[offBestIndex];
			}
			//if (!opVar.empty()) {
			//	if (opVar.back() - m_center.objective()[0] < -0.01) {
			//		std::cout << "center\t" << m_center.objective()[0] << std::endl;
			//		int stop = -1;
			//	}
			//}
			//opVar.push_back(m_center.objective()[0]);
		}
		for (auto& indi : offspring) {
			indi.setLocalSearch(false);
		}
		num_improve = 0;
		bool betterFlag(false);
		for (int pidx(0); pidx < pop.size(); ++pidx) {
			if (offspring[pidx].same(*pop[pidx], pro))continue;
			betterFlag = op.better(pro, offspring[pidx], *pop[pidx]);
			//rf!= op.evaluate(pro, alg, alg.idRandom(), offspring[pidx]);
			Real rPop(pop[pidx]->variableDistance(m_center, pro));
			Real rOff(offspring[pidx].variableDistance(m_center, pro));
			pop[pidx]->setImproved(false);
			if (rPop <= m_radius && rOff <= m_radius) {
				if(op.better(pro, offspring[pidx], *pop[pidx]))
				 {
					*pop[pidx] = offspring[pidx];
					pop[pidx]->setImproved(true);
					num_improve++;
				}
			}
			else {
			//	std::cout << "rPopRadius\t" << rPop << "\trOff\t" << rOff << std::endl;
				if (rOff < rPop) {
					*pop[pidx] = offspring[pidx];
					if (op.better(pro, offspring[pidx], *pop[pidx])) {
						pop[pidx]->setImproved(true);
						num_improve++;
					}
				}
			}
		}
		if (m_localSearch) {
			int localSearchId(-1);
			for (int pidx(0); pidx < pop.size(); ++pidx) {
				if (!pop[pidx]->flagLocalSearch()) {
					if (localSearchId == -1 || 
						op.better(pro, *pop[pidx], *pop[localSearchId])) {
						localSearchId = pidx;
					}
				}
			}
			if (localSearchId != -1) {
				std::vector<Real> origin_objs(pop[localSearchId]->objective());
				op.localSearch(*pop[localSearchId], alg.idRandom(), pro, alg, m_localSearchTimes, 0);
				if (DNStrategy.judgeChange(pro, origin_objs, pop[localSearchId]->objective())) {
					pop[localSearchId]->setLocalSearch(false);
				}
				else pop[localSearchId]->setLocalSearch(true);
			}

		}

		return rf;
	}


}

#endif
