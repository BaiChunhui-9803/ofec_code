
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

#ifndef OFEC_GL_ADAPTOR_SEQ_EVAL_H
#define OFEC_GL_ADAPTOR_SEQ_EVAL_H

#include"../../../template/framework/gl/gl_adaptor.h"
#include"../../../template/combination/sequence/sequence_algorithm.h"
#include "../../../../../utility/functional.h"
#include "../algo_util/fitness_calculator.h"
#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO

#include "../algo_util/density_matrix.h"

namespace ofec {

	template<typename TSequenceOperator>
	class AdaptorGLSeqEval: public AdaptorGL<typename TSequenceOperator::SolutionType> {

	public:
		//struct dataInfo {
		//	int m_id = 0;
		//	double m_dis2opt = 0;
		//	double m_dis2pop = 0;
		//	// map to [0,1]
		//	double m_fitness = 0;
		//};


		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;

	protected:
		std::vector<std::vector<Real>> m_pro;
		std::vector<std::vector<Real>> m_density_matrix;
		std::vector<double> m_fitness;
		FitnessCalculator m_fitness_calculator;
		DensityMatrix m_den_mat;
		
		

		enum class proUpdate { origin, min };
		proUpdate m_pro_update = proUpdate::min;
		int m_localSearchTimes = 30;
		Real m_radius = 1.0;
		SolutionType m_center;
		bool m_localSearch = true;


		//std::unique_ptr<EvaluationStrategy<SolutionType>> m_fitness_update;
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

#ifdef  OFEC_DEMO
		void updateDensityMatrix(
			const std::vector<std::vector<Real>> & pro,
			const std::vector<Real>& weight) {
			m_density_matrix = pro;
			double sum_weight(0);
			for (auto& it : weight) {
				sum_weight += it;
			}
			for (auto& it : m_density_matrix) {
				for (auto& it2 : it) {
					it2 /= sum_weight;
				}
			}
			
			//for (int idx(0); idx < m_density_matrix.size(); ++idx) {
			//	double sum(0);
			//	for (auto& nei : m_density_matrix[idx]) {
			//		sum += nei;
			//	}
			//	for (auto& nei : m_density_matrix[idx]) {
			//		nei /= sum;
			//	}
			//}

			//for (int idx(0); idx < m_density_matrix.size(); ++idx) {
			//	double maxValue(0);
			//	for (auto& nei : m_density_matrix[idx]) {
			//		maxValue =std::max(maxValue, nei);
			//	}
			//	if (maxValue > 0) {
			//		for (auto& nei : m_density_matrix[idx]) {
			//			nei /= maxValue;
			//		}
			//	}
			//}
		
		}
#endif
	public:
		AdaptorGLSeqEval() {
			m_alpha = 0.8;
		}
		virtual ~AdaptorGLSeqEval() = default;
		//double getPro(SolutionType& indi) {
		//	indi.calEdges();
		//}

		virtual void initialize(Problem *pro,Algorithm *alg) {
			m_den_mat.initialize(pro, alg);
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


		//virtual void updatePointDensity(Problem *pro, Algorithm *alg,
		//	std::vector<std::unique_ptr<SolutionType>>& pop) {
		////	using namespace ofec_demo;
		//	auto& optSol = pro->optBase().variable(0);
		//	std::vector<dataInfo> m_cur_sol_infos;
		//	auto& pointDensity = m_cur_sol_infos;
		//	pointDensity.resize(pop.size());
		//	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		//	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>);
		//	auto obj_range(pro->optBase().objRange());
		//	for (int idx(0); idx < pop.size(); ++idx) {
		//		//interpreter.updateEdges()
		//		auto& it(pointDensity[idx]);
		//		it.m_dis2opt = pop[idx]->variableDistance(optSol, pro);
		//		it.m_fitness = pop[idx]->objective()[0];
		//		if (it.m_fitness < obj_range.first) it.m_fitness = obj_range.first;
		//		else if (it.m_fitness > obj_range.second) it.m_fitness = obj_range.second;
		//		it.m_fitness = mapReal<double>(it.m_fitness, obj_range.first, obj_range.second, 0, 1);
		//		if (pro->optimizeMode()[0] == OptimizeMode::kMinimize) {
		//			it.m_fitness = 1.0 - it.m_fitness;
		//		}
		//		it.m_id = idx;
		//		it.m_dis2pop = 0;
		//		for (auto& edge : pop[idx]->edges()) {
		//			it.m_dis2pop += m_density_matrix[edge.first][edge.second];
		//		}
		//		it.m_dis2pop /= m_density_matrix.size();
		//    }
		//}

		virtual void updateDensityMatrix(std::vector<std::unique_ptr<SolutionType>>& pop) {
			m_den_mat.clear();
			std::vector<double> fit(pop.size());
			for (int idx(0); idx < pop.size(); ++idx) {
				fit[idx] = m_fitness_calculator.ObjectiveToFitness(*pop[idx]);
				fit[idx] = m_fitness_calculator.getPopFitness(fit[idx]);
			}
			m_den_mat.update(pop, fit);
		}


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
				interpreter.updateEdges(pro,sol);
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
			double dis(optSol.variableDistance(optSol,pro));
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
	inline void AdaptorGLSeqEval<TSequenceOperator>::shrinkPop(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop)
	{
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
		for (auto& it : pop) {
			if (op.learn_from_other(*it, rnd, pro, interpreter, m_center, m_radius)) {
				//	std::cout << "distance2center\t" << m_center.variableDistance(*it,pro) << std::endl;
				it->reset();
				interpreter.stepFinal(pro, *it);
				op.evaluate(pro, alg, rnd, *it);
			}
		}
		//updateCenter(pro, pop);
	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeqEval<TSequenceOperator>::updateProbability
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
			for (int pidx(0); pidx < weight.size(); ++pidx) {
				auto& indi(pop[pidx]);
				interpreter.updateEdges(pro, *indi);
				for (auto& edge : indi->edges()) {
					m_pro[edge.first][edge.second] += weight[pidx];
				}
			}
		}
	//	updateDensityMatrix(pop);
		if (m_pro_update == proUpdate::min) {
			modifyWeightMin(weight);
		}
		//auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());

	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeqEval<TSequenceOperator>::updateProbability(Problem *pro, Algorithm *alg, std::vector<SolutionType>& pop, const std::vector<Real>& weight, const std::vector<int>* index)
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
		//updateDensityMatrix(pop);

		if (m_pro_update == proUpdate::min) {
			modifyWeightMin(weight);
		}
		//auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());

	}


	template<typename TSequenceOperator>
	inline void AdaptorGLSeqEval<TSequenceOperator>::createSolution(Problem *pro, Algorithm *alg, Random *rnd, std::vector<std::unique_ptr<SolutionType>>& pop, std::vector<SolutionType>& offspring)
	{
		auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());

		for (int offId(0); offId < offspring.size(); ++offId)
		{
			auto& it = offspring[offId];
			it.reset();
			interpreter.stepInit(pro, it);
			while (!interpreter.stepFinish(pro, it)) {
				int next(-1);
				if (interpreter.stepFeasible(pro, it)) {
					if (rnd->uniform.next() < m_alpha
						/*std::max(m_alpha,1.0-m_radius/it->variable().size())*/) {
						next = op.learn_from_other(it, rnd, pro, interpreter, *pop[offId]);
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
	inline int AdaptorGLSeqEval<TSequenceOperator>::updateSolution(Problem *pro, Algorithm *alg, std::vector<std::unique_ptr<SolutionType>>& pop, std::vector<SolutionType>& offspring, int& num_improve)
	{

		auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
		int rf = 0;
		for (int pidx(0); pidx < pop.size(); ++pidx) {
			//	if (offspring[pidx].same(*pop[pidx], pro))continue;
			rf |= op.evaluate(pro, alg, alg.idRandom(), offspring[pidx]);
			rf |= op.evaluate(pro, alg, alg.idRandom(), *pop[pidx]);
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
				if (op.better(pro, offspring[pidx], *pop[pidx]))
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
				pop[localSearchId]->setLocalSearch(true);
			}

		}

		return rf;
	}
}

#endif



