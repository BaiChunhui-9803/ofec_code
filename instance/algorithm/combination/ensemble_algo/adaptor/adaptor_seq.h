/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yiya Diao
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
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
* Operator adaptator of genetic learning algorithm
*
*********************************************************************************/

#ifndef OFEC_ADAPTOR_SEQ_H
#define OFEC_ADAPTOR_SEQ_H

#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO

#include "../../../template/framework/gl/gl_adaptor.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"
#include "../../../template/combination/sequence/sequence_action.h"

#include "../../../combination/ensemble_algo/ensemble/sequence_actions_ensemble.h"

namespace ofec {

	template<typename TSequenceOperator>
	class AdaptorSeq {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;

	protected:
		std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_stra;
		int m_id_pop = 0;
		std::unique_ptr<SequenceAction<OpType>> m_step_action;
		int m_localSearchTimes = 30;
		bool m_localSearch = true;

		SolutionType m_center;
	public:
		AdaptorSeq() = default;
		virtual ~AdaptorSeq() = default;
		void setPopId(int popId) {
			m_id_pop = popId;
		}
		void setEvalStrategy(const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra) {
			m_eval_stra = eval_stra;
		}

		virtual void initialize(Problem *pro, Algorithm *alg);

		virtual void updateMemory( Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) {
			return m_step_action->globalUpdate(pro,rnd,  pop);
		}
		virtual void createSolution(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring);
		virtual int updateSolution(Problem *pro, Algorithm *alg,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring, int& num_improve);
#ifdef OFEC_DEMO

		virtual void updateOptPointInfo(Problem *pro, Algorithm *alg) {

			auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
			auto& optBase = pro->optBase();
			auto& optSol = pro->optBase().variable(0);
			auto& pointDensity = ofec_demo::BufferPointDistribution::m_opt_infos;
			pointDensity.resize(optBase.numberVariables());
			for (int idx(0); idx < optBase.numberVariables(); ++idx) {
				//interpreter.updateEdges()
				auto& it(pointDensity[idx]);
				SolutionType sol(optBase.variable(idx));
				sol.evaluate(pro, alg, false);
				interpreter.updateEdges(pro, sol);
				it.m_dis2opt = sol.variableDistance(optSol, pro);
				it.m_fitness = 1.0;
				it.m_id = idx;
				double fit = 1.0;
				it.m_dis2pop = m_den_mat.disToPop(sol, m_den_mat.getPopFitness(fit));
			}
		}

		virtual void updatePointDensity(Problem *pro, Algorithm *alg,
			std::vector<std::unique_ptr<SolutionType>>& pop) {
			using namespace ofec_demo;
			auto& optSol = pro->optBase().variable(0);
			double dis(optSol.variableDistance(optSol, pro));
			auto& pointDensity = ofec_demo::BufferPointDistribution::m_cur_sol_infos;
			pointDensity.resize(pop.size());
			//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
			//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>);
			auto obj_range(pro->optBase().objRange());
			for (int idx(0); idx < pop.size(); ++idx) {
				//interpreter.updateEdges()
				auto& it(pointDensity[idx]);
				it.m_dis2opt = pop[idx]->variableDistance(optSol, pro);
				it.m_fitness = pop[idx]->objective()[0];
				if (it.m_fitness < obj_range.first) it.m_fitness = obj_range.first;
				else if (it.m_fitness > obj_range.second) it.m_fitness = obj_range.second;
				it.m_fitness = mapReal<double>(it.m_fitness, obj_range.first, obj_range.second, 0, 1);
				if (pro->optimizeMode()[0] == OptimizeMode::kMinimize) {
					it.m_fitness = 1.0 - it.m_fitness;
				}
				it.m_id = idx;
				double fit = m_den_mat.ObjectiveToFitness(*pop[idx]);
				it.m_dis2pop = m_den_mat.disToPop(*pop[idx], m_den_mat.getPopFitness(fit));

				//it.m_dis2pop = 0;
				//for (auto& edge : pop[idx]->edges()) {
				//	it.m_dis2pop += m_density_matrix[edge.first][edge.second];
				//}
				//it.m_dis2pop /= pro->numberVariables();
				//it.m_dis2pop = 1.0 - it.m_dis2pop;
			}
		}
#endif
	};


	template<typename TSequenceOperator>
	inline void AdaptorSeq<TSequenceOperator>::initialize(Problem *pro, Algorithm *alg) {
		SequenceActionEnsembleParametersBase par;
		par.initGL_GA(pro, alg);
		m_step_action.reset(new SequenceActionEnsemble<TSequenceOperator>());
		m_step_action->initialize(par);
	}

	template<typename TSequenceOperator>
	inline void AdaptorSeq<TSequenceOperator>::createSolution(
		Problem *pro, Algorithm *alg, Random *rnd,
		std::vector<SolutionType*>& pop,
		std::vector<SolutionType>& offspring) {

		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());


		for (int offId(0); offId < offspring.size(); ++offId)
		{
			auto& it = offspring[offId];
			it.reset();
			interpreter.stepInit(pro, it);
			while (!interpreter.stepFinish(pro, it)) {
				int next = m_step_action->learn(it,rnd, pro, interpreter, op, pop, offId);
				if (next != -1) {
					interpreter.stepNext(pro, it, next);
				}
				else {
					interpreter.stepBack(pro, it);
				}
			}
			//	op.learn_from_other(it, rnd, pro, interpreter, m_center, m_radius);
			interpreter.stepFinal(pro, it);
		}
	}

	//template<typename TSequenceOperator>
	//inline void AdaptorGLSeqEval<TSequenceOperator>::createSolution(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop, std::vector<SolutionType>& offspring)
	//{
	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
	//	auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());

	//	for (int offId(0); offId < offspring.size(); ++offId)
	//	{
	//		auto& it = offspring[offId];
	//		it.reset();
	//		interpreter.stepInit(pro, it);
	//		while (!interpreter.stepFinish(pro, it)) {
	//			int next(-1);
	//			if (interpreter.stepFeasible(pro, it)) {
	//				if (rnd->uniform.next() < m_alpha
	//					/*std::max(m_alpha,1.0-m_radius/it->variable().size())*/) {
	//					next = op.learn_from_other(it, rnd, pro, interpreter, *pop[offId]);
	//				}
	//				if (next == -1) {
	//					auto weightFun = [&](const SolutionType& indi, int to) {
	//						return m_pro[interpreter.curPositionIdx(pro, indi)][to];
	//					};
	//					next = op.learn_from_global(it, rnd, pro, interpreter, weightFun);
	//				}
	//			}
	//			if (next != -1) {
	//				interpreter.stepNext(pro, it, next);
	//			}
	//			else {
	//				interpreter.stepBack(pro, it);
	//			}
	//		}
	//		op.learn_from_other(it, rnd, pro, interpreter, m_center, m_radius);
	//		interpreter.stepFinal(pro, it);
	//	}

	//}
	template<typename TSequenceOperator>
	inline int AdaptorSeq<TSequenceOperator>::updateSolution(
		Problem *pro, Algorithm *alg, std::vector<SolutionType*>& pop, std::vector<SolutionType>& offspring, int& num_improve)
	{

		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
		int rf = 0;
		// getBestIndex
		{

			int popBestIndex = calIdx<SolutionType*>(
				pop, [&]( SolutionType* a, SolutionType* b) {
				return  a->fitness() > b->fitness();
			});
			int offBestIndex = calIdx<SolutionType>(offspring,
				[&](const SolutionType& a, const SolutionType& b) {
				return a.fitness() > b.fitness();
			});


			if (pop[popBestIndex]->fitness() > offspring[offBestIndex].fitness()) {
				m_center = *pop[popBestIndex];
			}
			else {
				m_center = offspring[offBestIndex];
			}
		}
		for (auto& indi : offspring) {
			indi.setLocalSearch(false);
		}
		num_improve = 0;
		bool betterFlag(false);
		for (int pidx(0); pidx < pop.size(); ++pidx) {
			if (offspring[pidx].same(*pop[pidx], pro))continue;
			betterFlag = offspring[pidx].fitness() > pop[pidx]->fitness();
			Real rPop(pop[pidx]->variableDistance(m_center, pro));
			Real rOff(offspring[pidx].variableDistance(m_center, pro));
			pop[pidx]->setImproved(false);
			if (rOff < rPop) {
				*pop[pidx] = offspring[pidx];
				if (offspring[pidx].fitness() > pop[pidx]->fitness()) {
					pop[pidx]->setImproved(true);
					num_improve++;
				}
			}
		}
		if (m_localSearch) {
			int localSearchId(-1);
			for (int pidx(0); pidx < pop.size(); ++pidx) {
				if (!pop[pidx]->flagLocalSearch()) {
					if (localSearchId == -1 ||
						pop[pidx]->fitness() > pop[localSearchId]->fitness()) {
							{
								localSearchId = pidx;
							}
					}
				}
				if (localSearchId != -1) {

				
					rf |= op.localSearch(*pop[localSearchId], m_localSearchTimes,
						pro, alg, alg.idRandom(),
						m_eval_stra, m_id_pop
					);
					pop[localSearchId]->setLocalSearch(true);
				}
			}

			return rf;
		}
	}
}
#endif