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

#ifndef OFEC_POP_UNCERTIANTY_SEQ_ENSEMBLE_H
#define OFEC_POP_UNCERTIANTY_SEQ_ENSEMBLE_H


#ifdef OFEC_DEMO
#include <core/global_ui.h>
#include <ui/custom/buffer/algorithm/combination/multi_point_distribution/buffer_multi_point_distribution.h>
#endif
#include "../ensemble_algo/ensemble/sequence_actions_ensemble.h"
#include "../../template/framework/uncertianty/pop_uncertianty.h"
#include "../../../../utility/functional.h"



namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class PopulationUncertiantySeqEnsemble :
		virtual public PopUncertianty<TDistanceCalculator, typename TSequenceOperator::SolutionType> {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename TSequenceOperator::SolutionType;
			using InterpreterType = typename TSequenceOperator::InterpreterType;

		public:

			virtual void setId(int id)override {
				PopUncertianty::setId(id);
				m_actions.setPopId(id);
			}
			virtual void setEvalStrategy(
				const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra)override {
				PopUncertianty::setEvalStrategy(eval_stra);
				m_actions.setEvalStrategy(eval_stra);
			}

			virtual void initParamenters(Problem *pro, Algorithm *alg)override {
				PopUncertianty::initParamenters(pro, alg);
				SequenceActionEnsembleParametersBase par;
				par.initGL_GA_LS(pro, alg);
				m_actions.initialize(par);
				m_updateIter = 0;
				m_cur_iters = 0;
				m_converge_radius = 1.0;
				m_last_best.reset(nullptr);
				m_representative_id = -1;
				
			}

			virtual void updateIndisInfo() {
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				for (auto& it: m_individuals) {
					interpreter.updateEdges(m_problem.get() ,*it);
				}
				for (auto& it : m_offsprings) {
					interpreter.updateEdges(m_problem.get(), it);
				}
			}

			virtual void updateMemory()override;

			virtual void generateOffstrpings()override;

			//virtual int selectSurvivors() override;


			virtual void initialize(const std::vector<std::unique_ptr<SolutionType>>& inds) {
				resize(inds.size(), m_problem.get(), m_problem->numberVariables());
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				auto& op(GET_ASeq(this).getOp<OpType>());
				int idx(0);
				for (auto it(inds.begin()); it != inds.end(); ++it) {
					setSolution(idx, **it);
					m_individuals[idx]->reset();
					interpreter.stepFinal(m_problem.get(), *m_individuals[idx]);
					interpreter.updateEdges(m_problem.get(), *m_individuals[idx]);
					++idx;
				}
				updateSurvivors();
				m_distance_calculator->updateMemory(m_survivors);
				m_eval_strategy->calFitness(m_survivors, m_id);
				m_offsprings.resize(m_individuals.size());
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_offsprings[idx] = *m_individuals[idx];
				}
				//resize(temp_inds.size(), m_problem.get());

				resetCurState();
			}


			virtual void initialize(const std::list<std::unique_ptr<SolutionType>>& inds) {
				resize(inds.size(), m_problem.get(),m_problem->numberVariables());
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				auto& op(GET_ASeq(this).getOp<OpType>());
				int idx(0);
				for (auto it(inds.begin()); it != inds.end(); ++it) {
					setSolution(idx, **it);
					m_individuals[idx]->reset();
					interpreter.stepFinal(m_problem.get(), *m_individuals[idx]);
					interpreter.updateEdges(m_problem.get(), *m_individuals[idx]);
					++idx;
				}
				updateSurvivors();
				m_distance_calculator->updateMemory(m_survivors);
				m_eval_strategy->calFitness(m_survivors, m_id);
				m_offsprings.resize(m_individuals.size());
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_offsprings[idx] = *m_individuals[idx];
				}


				resetCurState();
				//resize(temp_inds.size(), m_problem.get());

			}

			virtual void reduceByFitness(int maxPop) {
				auto& pops(m_individuals);
				if (pops.size() <= m_init_popsize)return;
				std::vector<int> fitness_order;
				std::vector<int> belong_idx;
				std::vector<int> dis_order;
				sortByFitnessDis(pops, fitness_order, belong_idx, dis_order, m_problem.get());
				for (int idx(pops.size() - 1); idx >= maxPop; --idx) {
					int cur_id(dis_order[idx]);
					m_individuals[cur_id].reset(nullptr);
				}
				int cur_size(0);
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					if (m_individuals[idx] != nullptr) {
						if (idx != cur_size)
							m_individuals[cur_size++] = std::move(m_individuals[idx]);
					}
				}
				resize(cur_size, m_problem.get(), m_problem->numberVariables());
			}

			virtual void initialize(Problem *pro, Random *rnd) {

				PopUncertianty::initialize(pro, m_random.get());
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				auto& op(GET_ASeq(this).getOp<OpType>());
				auto& inds(m_individuals);
				int idx(0);
				for (auto it(inds.begin()); it != inds.end(); ++it) {
					//setSolution(idx, **it);
					m_individuals[idx]->reset();
					interpreter.stepFinal(m_problem.get(), *m_individuals[idx]);
					interpreter.updateEdges(m_problem.get(), *m_individuals[idx]);
					++idx;
				}
				updateSurvivors();
				m_eval_strategy->calFitness(m_survivors, m_id);
				m_offsprings.resize(m_individuals.size());
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_offsprings[idx] = *m_individuals[idx];
				}

				resetCurState();
			}


			template<typename ... Args>
			void resize(size_t size, Problem *pro, Args&& ...args){
				int origin_size(m_individuals.size());
				PopUncertianty::resize(size, pro, std::forward<Args>(args)...);
				auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
				if (size > origin_size) {
					for (int idx(origin_size); idx < m_individuals.size(); ++idx) {
						auto& it(m_individuals[idx]);
						it->reset();
						interpreter.stepFinal(pro, *it);
						interpreter.updateEdges(pro, *it);
					}
				}
			}


			virtual int insertNewIndis(int numPops) {
				resize(numPops,m_problem.get());
				return m_eval_strategy->calFitness(m_survivors, m_id);
			}

		
			virtual void deleteIndisFit(int numPop) {
				int originSize(m_individuals.size());
				std::sort(m_individuals.begin(), m_individuals.end(), [](
					const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b
					) {
					return a->fitness() > b->fitness();
				});

				m_individuals.resize(numPop);
				for (int idx(0); idx < originSize; ++idx) {
					m_individuals[idx]->setId(idx);
				}
			}
			void updateRrepresentative();

			virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
				int rf = 0;
				rf = PopUncertianty::evolve(pro, alg, rnd);
	
				updateConvergeState();
				return rf;
			}

			virtual void updateConvergeState() {
				m_representative_id = -1;
				getRepresentative();
				updateSurvivors();
				m_distance_calculator->updateMemory(m_survivors);
				++m_cur_iters;
				if (m_last_best == nullptr) {
					updateCurrentState();
				}
				else if (!m_last_best->same(getRepresentative(),m_problem.get()) || m_converge_radius > m_distance_calculator->avgRadius()) {
					updateCurrentState();
				}	
			}

			virtual void updateCurrentState() {
				m_updateIter = m_cur_iters;
				m_converge_radius = m_distance_calculator->avgRadius();
				m_last_best.reset(new SolutionType(getRepresentative()));

			}

			virtual void resetCurState() {
				updateSurvivors();
				updateRrepresentative();
				m_distance_calculator->updateMemory(m_survivors);
				updateCurrentState();
			}

			const SolutionType& getRepresentative() {
				if (m_representative_id == -1) {
					updateRrepresentative();
				}
				return *m_individuals[m_representative_id];
			}

			const SolutionType& getRepresentative() const {
				int representative_id = 0;
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					if (m_individuals[representative_id]->fitness() < m_individuals[idx]->fitness()) {
						representative_id = idx;
					}
				}
				return *m_individuals[representative_id];
			}


			virtual bool judgeConverged() const {
				return m_cur_iters - m_updateIter >= m_NotConvergeTimeThreadhold;
			}

			virtual void setInitRadius(Real radius) {
				m_initRadius = radius;
			}
			virtual Real getInitRadius()const {
				return m_initRadius;
			}
			
		protected:
			int m_representative_id = -1;
			int m_NotConvergeTimeThreadhold = 50;
			int m_updateIter = 0;
			int m_cur_iters = 0;
			double m_converge_radius = 1.0;
			std::unique_ptr<SolutionType> m_last_best=nullptr;


			SequenceActionEnsemble<OpType> m_actions;

			Real m_initRadius = 0.0;
			int m_init_popsize = 50;

	    public:

#ifdef  OFEC_DEMO
			static std::list<SolutionType> m_his_opts;
			static int m_curTime;
			static const int m_maxOpt;

			static void initHisInfo() {
				m_his_opts.clear();
				m_curTime = -1;
			}

			static void updateHisOpts(Problem *pro,Algorithm *alg) {
				
				
				auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
				//auto& interpreter = getInterpreter<InterpreterType>();
				auto& optSol = pro->optBase().variable(0);
				int curTime(GET_DOP(pro)->getCurTime());
				if (m_curTime != curTime) {
					m_curTime = curTime;
					m_his_opts.push_front(SolutionType(optSol));
					interpreter.updateEdges(pro, m_his_opts.front());
					if (m_his_opts.size() > m_maxOpt) {
						m_his_opts.pop_back();
					}
				}

				auto& opt_offset(ofec_demo::BufferMultiPointDistribution::m_opt_offset);
				if (m_his_opts.size() >= 2) {
					auto cur(m_his_opts.begin());
					auto before(m_his_opts.begin());
					++before;
					opt_offset = cur->variableDistance(*before, pro);
				}
				else {
					opt_offset = 1.0;
				}
				{
					int curId(m_curTime);
					for (auto& it : m_his_opts) {
						it.evaluate(pro, alg, false);
						it.setId(curId--);
					}
				}
			}


			static void updateDisInfo(
				ofec_demo::BufferMultiPointDistribution::SolDataInfo& solInfo,
				const PopulationUncertiantySeqEnsemble* pop, const SolutionType& sol,
				const std::pair<Real, Real>& objRange, SolutionType& opt_sol,
				Problem *pro
			);


			static void updateBestSolInfo(
				std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo>& solInfo,
				std::vector<const SolutionType*> sols,
				const std::pair<Real, Real>& obj_range,
				Problem *pro
			);

			static void updateBestSolInfo(
				std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo*>& solInfo,
				std::vector<const SolutionType*> sols,
				const std::pair<Real, Real>& obj_range,
				Problem *pro
			);


			static void setPopSolsInfo(
				std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo*>& solInfo,
				int from_sol_id,
				const PopulationUncertiantySeqEnsemble* pop,
				int numPop,
				Problem *pro
			);


			static void getFitnessOrder(
				const PopulationUncertiantySeqEnsemble* pop,
				std::vector<int>& fit_order
			);


			static void drawMultiPop(
				const std::vector< const PopulationUncertiantySeqEnsemble*> &sub_pop,
				Problem *pro,int show_offset=1
			);


			//static void updateGlobal


#endif

	};
#ifdef  OFEC_DEMO

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	std::list<typename TSequenceOperator::SolutionType> PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::m_his_opts;
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	int PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::m_curTime = -1;
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	const int PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::m_maxOpt = 5;


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::updateDisInfo(
			ofec_demo::BufferMultiPointDistribution::SolDataInfo& cur_sol_info,
			const PopulationUncertiantySeqEnsemble* pop, const SolutionType& sol,
			const std::pair<Real, Real>& obj_range, SolutionType& opt_sol,
			Problem *pro
			) {
		cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->disToPop(sol);
		cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->normalize01(cur_sol_info.m_dis2pop);

		cur_sol_info.m_dis2opt = opt_sol.variableDistance(sol, pro);
		cur_sol_info.m_objectives = mapReal<double>
			(sol.objective()[0], obj_range.first, obj_range.second,
				0, 1);
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::updateBestSolInfo(
			std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo>& best_sol,
			std::vector<const SolutionType*> best_indis,
			const std::pair<Real, Real>& obj_range,
			Problem *pro
		) {
		std::vector<int> fitness_order;
		std::vector<int> belong_idx;
		std::vector<int> distance_order;
		sortByFitnessDis(best_indis, fitness_order, belong_idx, distance_order, pro);

		for (int idx(0); idx < best_indis.size(); ++idx) {
			best_sol[idx].m_objectives = best_indis[idx]->objective()[0];
		//	if(best_sol[idx].m_objectives> obj_range.second) 
			best_sol[idx].m_objectives = 1.0 - mapRealInside<double>
				(best_sol[idx].m_objectives, obj_range.first, obj_range.second,
					0, 1);
			if (belong_idx[idx] == idx) {
				best_sol[idx].m_min_dis = 1.0;
			}
			else {
				best_sol[idx].m_min_dis = best_indis[idx]->variableDistance(*best_indis[belong_idx[idx]], pro);
			}
		}
	}



	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::updateBestSolInfo(
			std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo*>& best_sol,
			std::vector<const SolutionType*> best_indis,
			const std::pair<Real, Real>& obj_range,
			Problem *pro
		) {

		if (best_indis.empty())return;
		std::vector<int> fitness_order;
		std::vector<int> belong_idx;
		std::vector<int> distance_order;
		sortByFitnessDis(best_indis, fitness_order, belong_idx, distance_order, pro);



		for (int idx(0); idx < best_indis.size(); ++idx) {
			best_sol[idx]->m_objectives = best_indis[idx]->objective()[0];
			best_sol[idx]->m_objectives = 1.0 - mapRealInside<double>
				(best_sol[idx]->m_objectives, obj_range.first, obj_range.second,
					0, 1);
			if (belong_idx[idx] == idx) {
				best_sol[idx]->m_min_dis = 1.0;
			}
			else {
				best_sol[idx]->m_min_dis = best_indis[idx]->variableDistance(*best_indis[belong_idx[idx]], pro);
			}
		}
	
		double id(0);
		double total_id(fitness_order.size()-1);
		for (auto& it : fitness_order) {
			best_sol[it]->m_fitness_order_ratio = id / total_id;
			id = id + 1;
		}
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		setPopSolsInfo(
			std::vector<ofec_demo::BufferMultiPointDistribution::SolDataInfo*>& solInfo,
			int from_sol_id,
			const PopulationUncertiantySeqEnsemble* pop,
			int numPop,
			Problem *pro
		) {
		for (int bidx(0); bidx < pop->size(); ++bidx) {
			auto cur_sol(solInfo[bidx + from_sol_id]);
			cur_sol->m_show_pop_id = 0;
			cur_sol->m_pop_id = pop->id();
			cur_sol->m_sol_id = bidx;
			cur_sol->m_pop_ratio = double(cur_sol->m_pop_id) / double(numPop-1);
			cur_sol->m_sol_id_ratio = double(cur_sol->m_sol_id) / double(pop->size() - 1);
		}
	}


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		getFitnessOrder(
			const PopulationUncertiantySeqEnsemble* pop,
			std::vector<int>& fit_order
		) {
			{
				std::vector<int> order;
				order.resize(pop->size());
				for (int idx(0); idx < pop->size(); ++idx) {
					order[idx] = idx;
				}
				std::sort(order.begin(), order.end(),
					[&](int a, int b) {
					return (*pop)[a].fitness() > (*pop)[b].fitness();
				});

				fit_order.resize(order.size());
				for (int idx(0); idx < fit_order.size(); ++idx) {
					fit_order[order[idx]] = idx;
				}
			}
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		drawMultiPop(
			const std::vector< const PopulationUncertiantySeqEnsemble*>& sub_pop,
			Problem *pro, int show_offset
		) {

		auto& pop_idxs(ofec_demo::BufferMultiPointDistribution::m_popIdxs);
		pop_idxs.resize(sub_pop.size() + show_offset);
		pop_idxs[0] = -1;
		for (int idx(0); idx < sub_pop.size(); ++idx) {
			pop_idxs[idx+ show_offset] = sub_pop[idx]->id();
		}

		auto& popDivLines(ofec_demo::BufferMultiPointDistribution::m_pop_div_val);
		popDivLines.resize(sub_pop.size() + show_offset);
		popDivLines[0] = -1;
		for (int idx(0); idx < sub_pop.size(); ++idx) {
			popDivLines[idx+ show_offset] = sub_pop[idx]->getDistanceCalculator()->normalize01(sub_pop[idx]->getDistanceCalculator()->radiusThreadhold());
		}


		auto& obj_range(pro->optBase().objRange());
			{
				using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;
				SolDataInfo cur;

				auto& solInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_sols;
				auto& obestInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_obest;
				auto& optInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_opt;

				solInfo.clear();
				obestInfo.clear();
				optInfo.clear();

				const PopulationUncertiantySeqEnsemble* pop = nullptr;
				auto& opt_sol(m_his_opts.front());
				std::vector<int> fit_order;

				for (int idpop(0); idpop < sub_pop.size(); ++idpop) {
					cur.m_show_pop_id = idpop + show_offset;
					cur.m_opt_id = 0;
					cur.m_pop_id = idpop + 1;
					cur.m_pop_ratio = double(cur.m_pop_id) / double(sub_pop.size());
					pop = sub_pop[idpop];
					{
						std::vector<int> order;
						order.resize(pop->size());
						for (int idx(0); idx < pop->size(); ++idx) {
							order[idx] = idx;
						}
						std::sort(order.begin(), order.end(),
							[&](int a, int b) {
							return (*pop)[a].fitness() > (*pop)[b].fitness();
						});

						fit_order.resize(order.size());
						for (int idx(0); idx < fit_order.size(); ++idx) {
							fit_order[order[idx]] = idx;
						}
					}

					for (int solid(0); solid < pop->size(); ++solid) {
						cur.m_sol_id = solid;
						cur.m_sol_id_ratio = double(solid) / pop->size();
						cur.m_fitness_order_ratio = double(fit_order[solid]) / pop->size();
						auto& sol((*pop)[solid]);
						updateDisInfo(cur, pop, sol, obj_range, opt_sol, pro);
						solInfo.push_back(cur);
					}

					for (int popidy(0); popidy < sub_pop.size(); ++popidy) {
						if (popidy != idpop) {
							cur.m_pop_id = popidy + 1;
							cur.m_pop_ratio = double(cur.m_pop_id) / sub_pop.size();
							cur.m_sol_id = 0;
							cur.m_sol_id_ratio = 0;
							cur.m_fitness_order_ratio = 0;
							auto& sol(sub_pop[popidy]->getRepresentative());
								updateDisInfo(cur, pop, sol, obj_range, opt_sol, pro);
							obestInfo.push_back(cur);
						}
					}
					cur.m_pop_id = -1;
					cur.m_fitness_order_ratio = 0;
					for (auto solIter(m_his_opts.begin()); solIter != m_his_opts.end(); ++solIter) {
						cur.m_sol_id = solIter->id();
						cur.m_sol_id_ratio = double(cur.m_sol_id) / double(m_his_opts.size());
						auto& sol(*solIter);
						updateDisInfo(cur, pop, sol, obj_range, opt_sol, pro);
						optInfo.push_back(cur);
					}
				}
			}

	}


#endif



	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::updateMemory() {
		++m_memory_iter;
		updateSurvivors();
		
		auto& interpreter(GET_ASeq(this).getInterpreter<InterpreterType>());
		for (auto& it : m_survivors) {
			interpreter.updateEdges(m_problem.get(), *it);
		}

		m_distance_calculator->updateMemory(m_survivors);
		m_eval_strategy->updateMemory(m_survivors, m_id);
		m_actions.globalUpdate(m_problem.get(), m_random.get(), m_survivors);

		updateRrepresentative();
	}


	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::generateOffstrpings() {
		m_actions.createSolution(m_problem.get(), this, m_random.get(),
			m_survivors, m_offsprings);
	}
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::updateRrepresentative()
	{
		m_representative_id = 0;
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			if (m_individuals[m_representative_id]->fitness() < m_individuals[idx]->fitness()) {
				m_representative_id = idx;
			}
		}
	}

}
#endif

