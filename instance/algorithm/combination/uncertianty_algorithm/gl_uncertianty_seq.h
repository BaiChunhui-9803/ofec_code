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


#ifndef GL_UNCERTIANTY_SEQ_H
#define GL_UNCERTIANTY_SEQ_H


#include "../../../../core/problem/uncertainty/dynamic.h"
#include "../../../../instance/problem/combination/selection_problem/selection_problem.h"

#include"../../template/combination/sequence/sequence_algorithm.h"
#include"../../template/combination/sequence/sequence_interpreter.h"
#include"../../template/combination/sequence/sequence_operator.h"
#include"../../template/combination/multi_population/distance_calculator.h"
#include"../../combination/ensemble_algo/gl_actions/weight_calculator.h"


#include "pop_uncertainty_seq_ensemble.h"
#include "pop_uncertianty_seq.h"


#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO


namespace ofec {
	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	class GeneticLearningUncertiantySeq : public AlgorithmSeq {
	public:
		using OpType = typename TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
	protected:
		PopulationUncertiantySeqEnsemble<TDistanceCalculator, OpType> m_pop;
		std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_strategy;
		virtual void initSolutions(int numIndis);
		
#ifdef  OFEC_DEMO
		void updateBuffer();
	//	void updateOptPointInfo();
	//	void updatePointDensity();
#endif //  OFEC_DEMO


	public:
		virtual void initialize_()override;
		virtual void run_()override;
		// Í¨¹ý Algorithm ¼Ì³Ð
		virtual void record() override {}


	protected:
		
#ifdef  OFEC_DEMO
		std::list<SolutionType> m_his_opts;
		int m_curTime = -1;
		const int m_maxOpt = 5;
#endif

	};
	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void GeneticLearningUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TSequenceOperator>::initSolutions(int numIndis) {

		m_pop.resize(numIndis, m_problem.get(), m_problem->numberVariables());
		m_pop.initialize(m_problem.get(), m_random.get());


	}

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void GeneticLearningUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TSequenceOperator>::initialize_()
	{
		AlgorithmSeq::initialize_();
		m_eval_strategy.reset(new  TEvaluationStrategy<SolutionType>());
		m_eval_strategy->initialize(m_problem.get(), this);
		m_keep_candidates_updated = true;
		assignInterpreter<InterpreterType>();
		assignOps<OpType>();
		m_pop.initParamenters(m_problem.get(), this);
		m_pop.setEvalStrategy(m_eval_strategy);
		m_pop.setId(0);
		auto& v = *m_param;
		initSolutions(std::get<int>(m_param->at("population size")));
		m_pop.updateMemory();
		//m_pop.assign(std::get<int>(m_param->at("population size")), m_problem.get());
		//m_pop.resize(std::get<int>(m_param->at("population size")), m_problem.get());
		//m_pop.initialize(m_problem.get(), this, m_random.get());


#ifdef  OFEC_DEMO
		m_his_opts.clear();
		m_curTime = -1;
#endif

	}

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void GeneticLearningUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TSequenceOperator>::run_()
	{
		while (!terminating()) {
			//	std::cout << "alg id\t" << this << std::endl;
			m_pop.evolve(m_problem.get(), this, m_random.get());
			m_pop.updateMemory();
#ifdef  OFEC_DEMO
			updateBuffer();
#endif

			//std::cout << m_pop.best().front()->objective()[0]- m_problem->optBase().objective(0)[0]<< std::endl;
			//m_pop.updateBest(m_problem.get());
		//	std::cout << evaluations() << std::endl;
		//	m_problem->showInfomations(this);
		}
	}

#ifdef  OFEC_DEMO
	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void GeneticLearningUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TSequenceOperator>::updateBuffer() {


		ofec_demo::BufferMultiPointDistribution::m_numOpts
			= m_problem->optBase().numberVariables();
		ofec_demo::BufferMultiPointDistribution::m_numPops
			= 1;
		int numPops(1);

		using sol_info_type = ofec_demo::BufferMultiPointDistribution::SolDataInfo;

		auto& pop_idxs(ofec_demo::BufferMultiPointDistribution::m_popIdxs);
		pop_idxs.resize(numPops);
		pop_idxs[0] = 0;

		auto& popDivLines(ofec_demo::BufferMultiPointDistribution::m_pop_div_val);
		popDivLines.resize(1);
		popDivLines[0] = m_pop.getDistanceCalculator()->normalize01(m_pop.getDistanceCalculator()->radiusThreadhold());

		auto& interpreter = getInterpreter<InterpreterType>();
		auto& optSol = m_problem->optBase().variable(0);
		int curTime(GET_DOP(m_problem.get())->getCurTime());
		if (m_curTime != curTime) {
			m_curTime = curTime;
			m_his_opts.push_front(SolutionType(optSol));
			interpreter.updateEdges(m_problem.get(), m_his_opts.front());
			if (m_his_opts.size() > m_maxOpt) {
				m_his_opts.pop_back();
			}
		}


		{
			int curId(m_curTime);
			for (auto& it : m_his_opts) {
				it.evaluate(m_problem.get(), this, false);
				it.setId(curId--);
			}
		}


		auto& obj_range(m_problem->optBase().objRange());
		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_cur_sols);
			sol_info.clear();
			sol_info_type cur_sol_info;
			auto& opt_sol(m_his_opts.front());
			const PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>* pop = &m_pop;
			{
				cur_sol_info.m_show_pop_id = 0;
				cur_sol_info.m_opt_id = 0;
				cur_sol_info.m_pop_id = 0;
				for (int solid(0); solid < pop->size(); ++solid) {
					cur_sol_info.m_sol_id = solid;
					auto& sol((*pop)[solid]);
					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());
					sol_info.push_back(cur_sol_info);
				}

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

					for (int idx(0); idx < order.size(); ++idx) {
						sol_info[order[idx]].m_fitness_order_ratio = double(idx) / double(order.size() - 1);
					}
				}
			}
		}


		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_other_pop_sols);
			sol_info.clear();
		}


		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_opt_infos);
			sol_info.clear();
			sol_info_type cur_sol_info;
			auto& opt_sol(m_his_opts.front());
			const PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>* pop = &m_pop;
			{
				cur_sol_info.m_show_pop_id = 0;
				cur_sol_info.m_opt_id = 0;
				cur_sol_info.m_pop_id = -1;
				int solid(0);
				for (auto solIter(m_his_opts.begin()); solIter != m_his_opts.end(); ++solIter) {
					cur_sol_info.m_sol_id = solid;
					auto& sol(*solIter);
					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());


					//cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->disToPop(sol);
					//cur_sol_info.m_dis2opt = opt_sol.variableDistance(sol, m_problem.get());
					//cur_sol_info.m_objectives = mapReal<double>
					//	(sol.objective()[0], obj_range.first, obj_range.second,
					//		0, 1);
					sol_info.push_back(cur_sol_info);
					solid++;
				}
			}
		}

		//{
		//	std::vector<const SolutionType*> best_indis(m_pop.getSurvivors().size());

		//	auto& best_sol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
		//	best_sol.resize(best_indis.size());
		//	{


		//		int bidx(0);
		//		for (; bidx < m_pop.getSurvivors().size(); ++bidx) {
		//			best_indis[bidx] = m_pop.getSurvivors()[bidx];
		//			best_sol[bidx].m_show_pop_id = -1;
		//			best_sol[bidx].m_pop_id = 0;
		//			best_sol[bidx].m_opt_id = 0;

		//		}



		//		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		//			updateBestSolInfo(best_sol, best_indis, obj_range, m_problem.get());
		//	}
		//}



		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_opt_infos);
			sol_info.clear();
			sol_info_type cur_sol_info;
			auto& opt_sol(m_his_opts.front());
			const PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>* pop = &m_pop;
			{
				cur_sol_info.m_show_pop_id = 0;
				cur_sol_info.m_opt_id = 0;
				cur_sol_info.m_pop_id = -1;
				int solid(0);
				for (auto solIter(m_his_opts.begin()); solIter != m_his_opts.end(); ++solIter) {
					cur_sol_info.m_sol_id = solid;
					auto& sol(*solIter);

					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());

					//cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->disToPop(sol);
					//cur_sol_info.m_dis2opt = opt_sol.variableDistance(sol, m_problem.get());
					//cur_sol_info.m_objectives = mapReal<double>
					//	(sol.objective()[0], obj_range.first, obj_range.second,
					//		0, 1);
					sol_info.push_back(cur_sol_info);
					solid++;
				}
			}
		}

		{
			using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;


			std::vector<SolDataInfo*> best_sol(m_pop.getSurvivors().size());
			std::vector<const SolutionType*> best_indis(best_sol.size());

			auto& maintsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
			auto& subsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_subpop);
			auto& optsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_opt);
			{
			int bidx(0);
				/*		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
							setPopSolsInfo(best_sol, 0, &m_main_pop, m_sub_pops.size() + 1, m_problem.get());
			*/			
			maintsol.resize(m_pop.getSurvivors().size());
			for (; bidx < m_pop.getSurvivors().size(); ++bidx) {
				best_indis[bidx] = m_pop.getSurvivors()[bidx];
				best_sol[bidx] = &maintsol[bidx];

				auto& sol_info(best_sol[bidx]);
				sol_info->m_show_pop_id = 0;
				sol_info->m_pop_id = 0;
				sol_info->m_sol_id = bidx;
				sol_info->m_pop_ratio = 0;
				sol_info->m_sol_id_ratio = double(sol_info->m_sol_id) / double(m_pop.size() - 1);


			}

			//subsol.resize(m_sub_pops.size());

			//for (int pidx(0); pidx < m_sub_pops.size(); ++pidx) {
			//	best_indis[bidx + pidx] = &m_sub_pops[pidx].getRepresentative();
			//	best_sol[bidx + pidx] = &subsol[pidx];

			//	auto& sol_info(best_sol[bidx + pidx]);
			//	sol_info->m_show_pop_id = 0;
			//	sol_info->m_pop_id = pidx;
			//	sol_info->m_sol_id = 0;
			//	sol_info->m_pop_ratio = double(sol_info->m_pop_id) / double(m_sub_pops.size());
			//	sol_info->m_sol_id_ratio = 0;
			//}

			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
				updateBestSolInfo(best_sol, best_indis, obj_range, m_problem.get());;

			optsol.resize(1);

			{
				auto& sol_info(optsol.back());
				sol_info.m_show_pop_id = 0;
				sol_info.m_pop_id = 0;
				sol_info.m_sol_id = 0;
				sol_info.m_pop_ratio = 0;
				sol_info.m_sol_id_ratio = 0;

				sol_info.m_objectives = 1.0;
				sol_info.m_min_dis = 1.0;
				sol_info.m_fitness_order_ratio = 0;
			}
			}
		}



		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

#endif
}

#endif 