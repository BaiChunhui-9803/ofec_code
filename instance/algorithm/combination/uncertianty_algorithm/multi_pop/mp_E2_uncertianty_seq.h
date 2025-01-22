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


#ifndef MP_E2_UNCERTIANTY_SEQ_H
#define MP_E2_UNCERTIANTY_SEQ_H




#include "pop_mp_main_uncertianty_seq.h"
#include "pop_mp_sub_uncertianty_seq.h"
#include "pop_mp_diversity_main_uncertianty_seq.h"
#include "mp_uncertianty_seq.h"
#include "../../../../../core/algorithm/multi_population.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"

#include<list>

#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/multi_point_distribution/buffer_multi_point_distribution.h>

#endif //  OFEC_DEMO


// for test
#include "../../../combination/sequence/SP/sp_operator.h"
#include "../../../template/framework/uncertianty/distance_calculator.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"



namespace ofec {


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class MultiPopE2UncertiantySeq : public AlgorithmSeq {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename OpType::SolutionType;
			using InterpreterType = typename OpType::InterpreterType;
			using MainPopType = typename PopulationMPDivMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>;
			using SubPopType = typename PopulationMPSubUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>;

		protected:

			MainPopType m_main_pop;
			MultiPopulation<SubPopType> m_sub_pops;
			std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_strategy;
			IdManager m_ids;


			int m_numHis = 1000;
			int m_numPop = 20;
			int m_numMainPop = 60;

			//bool m_flag_

#ifdef  OFEC_DEMO
			std::list<SolutionType> m_his_opts;
			int m_curTime = -1;
			
			const int m_maxOpt = 5;
#endif

		public:
#ifdef  OFEC_DEMO
			void updateBuffer();
			
			//void updateOptPointInfo();
			//void updatePointDensity();
#endif //  OFEC_DEMO


		public:
			virtual void initialize_()override;
			virtual void run_()override;
			// Í¨¹ý Algorithm ¼Ì³Ð
			virtual void record() override {}

		protected:
			virtual void appendPop(std::vector<std::unique_ptr<SolutionType>>& inds);
			virtual void addNewlyPop();
			virtual void filterPop();
	};

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::initialize_()
	{

#ifdef  OFEC_DEMO
		m_his_opts.clear();
		m_curTime = -1;

		
		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
			initHisInfo();
#endif

		m_ids.initialize();
		AlgorithmSeq::initialize_();
		m_eval_strategy.reset(new TEvaluationStrategy<SolutionType>());
		m_eval_strategy->initialize(m_problem.get(), this);
		m_keep_candidates_updated = true;
		assignInterpreter<InterpreterType>();
		assignOps<OpType>();
		m_main_pop.initParamenters(m_problem.get(), this);
		m_main_pop.setEvalStrategy(m_eval_strategy);
		m_main_pop.setId(m_ids.getId());
		m_main_pop.resize(m_numMainPop,
			m_problem.get(),
			m_problem->numberVariables());
		m_main_pop.initialize(m_problem.get(), m_random.get());
		m_main_pop.updateSurvivors();
		m_main_pop.updateMemory();
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::
		appendPop(std::vector<std::unique_ptr<SolutionType>>& inds)
	{

		std::unique_ptr<SubPopType> curPop(new SubPopType);
		curPop->initParamenters(m_problem.get(), this);
		curPop->setEvalStrategy(m_eval_strategy);
		curPop->initialize(inds);
		
		curPop->setId(m_ids.getId());
		m_eval_strategy->clearMemory(curPop->id());
		curPop->updateMemory();

		curPop->generateNewSolution();
		curPop->updateMemory();
		
		m_sub_pops.append(curPop);
	}

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::addNewlyPop()
	{

		if (m_main_pop.getHis().size() >= m_numHis) {

			auto& inds(m_main_pop.getHis());
			for (auto& indi : inds) {
				m_eval_strategy->evaluate(*indi, m_main_pop.id());
			}

			FDC<SolutionType> clu(inds, m_problem.get());
			clu.calFitDist();
			///clu.calculate_objective_distance();
			clu.clusterBySize(m_numPop);
			auto& cluster_re(clu.clusters());
			
			std::vector<std::unique_ptr<SolutionType>> new_indis;
			for (int cidx(0); cidx < cluster_re.size(); ++cidx) {
				new_indis.clear();
				for (auto& sid : cluster_re[cidx]) {
					new_indis.push_back(std::move(inds[sid]));
					inds[sid].reset(nullptr);
				}
				appendPop(new_indis);
			}

			m_main_pop.setSaveFlag(false);
		}
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::filterPop()
	{


		bool erase_flag(false);
		int popx_id(0);
		auto& pops(m_sub_pops.getPops());
		int cur_id(0);
		int relation_ship(-1);
		std::vector<int> belong(pops.size(), -1);
		//std::vector<double> fitness(pops.size(), -std::numeric_limits<double>::max());
		for (int idx(0); idx < pops.size(); ++idx) {
			for (int idy(0); idy < pops.size(); ++idy) {
				if (idx != idy) {
					relation_ship = pops[idx]->judgeRelationship(*pops[idy]);

					if (relation_ship < 0) {
						belong[idx] = idy;
						break;
					}
					else if (relation_ship > 0) {
						if (belong[idy] != -1) {
							if (pops[idx]->getRepresentative().fitness() >
								pops[belong[idy]]->getRepresentative().fitness()) {
								belong[idy] = idx;
							}
						}
						else {
							belong[idy] = idx;
						}
					}

				}

			}
		}

		for (int idx(0); idx < pops.size(); ++idx) {
			if (belong[idx] != -1 && pops[belong[idx]]->id() != -1) {
				pops[belong[idx]]->mergePopByFitnessDis(*pops[idx]);
				m_ids.insertId(pops[idx]->id());
				pops[idx]->setId(-1);
			}
		}
		for (auto it(m_sub_pops.begin()); it != m_sub_pops.end();) {
			if ((**it).id() == -1) {
				it = m_sub_pops.remove(it);
			}
			else ++it;
		}
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::run_()
	{

		int iter(0);
		while (!terminating()) {
			++iter;
			//	std::cout << "alg id\t" << this << std::endl;
			m_main_pop.evolve(m_problem.get(), this, m_random.get());
			m_main_pop.updateRrepresentative();

			addNewlyPop();

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].evolve(m_problem.get(), this, m_random.get());
			}
			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateRrepresentative();
			}

			filterPop();
			for (auto it(m_sub_pops.begin()); it != m_sub_pops.end(); ) {
				if ((**it).judgeConverged()) {
					m_main_pop.insert(**it));
					it = m_sub_pops.remove(it);
				}
				else ++it;
			}
			

			std::vector<std::unique_ptr<SolutionType>> indis;
			for (auto it(m_sub_pops.begin()); it != m_sub_pops.end(); ++it) {
				indis.emplace_back(new SolutionType((**it).getRepresentative()));	
			}
			m_main_pop.insert(indis);


			if (m_sub_pops.size() <=1&&!m_main_pop.getSaveFlag()) {
				m_main_pop.setSaveFlag(true);
			}


			m_main_pop.updateSurvivors();
			m_main_pop.updateMemory();
			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateSurvivors();
				m_sub_pops[idx].updateMemory();
			}

#ifdef  OFEC_DEMO
			updateBuffer();
			m_problem->updateLocalOptimaInfo();
#endif
			//std::cout << m_pop.best().front()->objective()[0] - m_problem->optBase().objective(0)[0] << std::endl;
			//m_pop.updateBest(m_problem.get());
		//	std::cout << evaluations() << std::endl;
		//	m_problem->showInfomations(this);
		}
	}

#ifdef  OFEC_DEMO
	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopE2UncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::updateBuffer() {

		m_main_pop.updateIndisInfo();
		for (int idx(0); idx < m_sub_pops.size(); ++idx) {
			m_sub_pops[idx].updateIndisInfo();
		}


		ofec_demo::BufferMultiPointDistribution::m_numOpts
			= m_problem->optBase().numberVariables();
		ofec_demo::BufferMultiPointDistribution::m_numPops
			= m_sub_pops.size() + 1;


		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
			updateHisOpts(m_problem.get(), this);

		int numPops(m_sub_pops.size());

		using sol_info_type = ofec_demo::BufferMultiPointDistribution::SolDataInfo;

		auto& pop_idxs(ofec_demo::BufferMultiPointDistribution::m_popIdxs);
		pop_idxs.resize(numPops + 1);
		pop_idxs[0] = -1;
		for (int idx(1); idx < pop_idxs.size(); ++idx) {
			pop_idxs[idx] = m_sub_pops[idx - 1].id();
		}

		auto& popDivLines(ofec_demo::BufferMultiPointDistribution::m_pop_div_val);
		popDivLines.resize(numPops + 1);
		popDivLines[0] = -1;
		for (int idx(1); idx < pop_idxs.size(); ++idx) {
			popDivLines[idx] = m_sub_pops[idx - 1].getDistanceCalculator()->normalize01(m_sub_pops[idx - 1].getDistanceCalculator()->radiusThreadhold());
		}

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

		auto& opt_offset(ofec_demo::BufferMultiPointDistribution::m_opt_offset);
		if (m_his_opts.size() >= 2) {
			auto cur(m_his_opts.begin());
			auto before(m_his_opts.begin());
			++before;
			opt_offset = cur->variableDistance(*before, m_problem.get());
		}
		else {
			opt_offset = 1.0;
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
				<TDistanceCalculator, TAdaptor>* pop = nullptr;


			for (int idx(0); idx < 1; ++idx) {
				if (idx == 0)  pop = &m_main_pop;
				else pop = &m_sub_pops[idx - 1];
				cur_sol_info.m_show_pop_id = idx;
				cur_sol_info.m_opt_id = 0;
				cur_sol_info.m_pop_id = idx;
				for (int solid(0); solid < pop->size(); ++solid) {
					cur_sol_info.m_sol_id = solid;
					auto& sol((*pop)[solid]);
					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
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
			using pop_type = PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>;
			std::vector<const pop_type* > sub_pops(m_sub_pops.size());

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				sub_pops[idx] = &m_sub_pops[idx];
			}
			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
				drawMultiPop(sub_pops, m_problem.get());


		}



		//{

		//	auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_other_pop_sols);
		//	sol_info.clear();
		//	sol_info_type cur_sol_info;
		//	auto& opt_sol(m_his_opts.front());
		//	const PopulationUncertiantySeqEnsemble
		//		<TDistanceCalculator, TAdaptor>* pop = nullptr;

		//	const PopulationUncertiantySeqEnsemble
		//		<TDistanceCalculator, TAdaptor>* pop2 = nullptr;
		//	for (int idx(1); idx < pop_idxs.size(); ++idx) {
		//		if (idx == 0)  pop = &m_main_pop;
		//		else pop = &m_sub_pops[idx - 1];
		//		cur_sol_info.m_show_pop_id = idx;
		//		cur_sol_info.m_opt_id = 0;



		//		for (int idy(1); idy < pop_idxs.size(); ++idy) {
		//			cur_sol_info.m_pop_id = idy;
		//			if (idx != idy) {
		//				if (idy == 0)  pop2 = &m_main_pop;
		//				else pop2 = &m_sub_pops[idy - 1];
		//				/*for (int solid(0); solid < pop2->size(); ++solid)*/ {
		//					cur_sol_info.m_sol_id = 0;
		//					auto& sol((*pop2).getRepresentative());
		//					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());
		//					//	PopulationUncertiantySeqEnsemble::updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());;


		//						//cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->disToPop(sol);
		//						////cur_sol_info.m_dis2pop=pop->getDistanceCalculator()->

		//						//cur_sol_info.m_dis2opt = opt_sol.variableDistance(sol, m_problem.get());
		//						//cur_sol_info.m_objectives = mapReal<double>
		//						//	(sol.objective()[0], obj_range.first, obj_range.second,
		//						//		0, 1);
		//					sol_info.push_back(cur_sol_info);
		//				}
		//			}
		//		}
		//	}
		//}


		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_opt_infos);
			sol_info.clear();
			sol_info_type cur_sol_info;
			auto& opt_sol(m_his_opts.front());
			const PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>* pop = &m_main_pop;
			{
				cur_sol_info.m_show_pop_id = 0;
				cur_sol_info.m_opt_id = 0;
				cur_sol_info.m_pop_id = -1;
				int solid(0);
				for (auto solIter(m_his_opts.begin()); solIter != m_his_opts.end(); ++solIter) {
					cur_sol_info.m_sol_id = solid;
					auto& sol(*solIter);
					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());
					sol_info.push_back(cur_sol_info);
					solid++;
				}
			}
		}

		{
			using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;


			std::vector<SolDataInfo*> best_sol(m_main_pop.getSurvivors().size() + m_sub_pops.size());
			std::vector<const SolutionType*> best_indis(best_sol.size());

			auto& maintsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
			auto& subsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_subpop);
			auto& optsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_opt);
			{
				int bidx(0);
				/*		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
							setPopSolsInfo(best_sol, 0, &m_main_pop, m_sub_pops.size() + 1, m_problem.get());
			*/			maintsol.resize(m_main_pop.getSurvivors().size());
			for (; bidx < m_main_pop.getSurvivors().size(); ++bidx) {
				best_indis[bidx] = m_main_pop.getSurvivors()[bidx];
				best_sol[bidx] = &maintsol[bidx];

				auto& sol_info(best_sol[bidx]);
				sol_info->m_show_pop_id = 0;
				sol_info->m_pop_id = m_main_pop.id();
				sol_info->m_sol_id = bidx;
				sol_info->m_pop_ratio = double(sol_info->m_pop_id) / double(m_sub_pops.size());
				sol_info->m_sol_id_ratio = double(sol_info->m_sol_id) / double(m_main_pop.size() - 1);


			}

			subsol.resize(m_sub_pops.size());

			for (int pidx(0); pidx < m_sub_pops.size(); ++pidx) {
				best_indis[bidx + pidx] = &m_sub_pops[pidx].getRepresentative();
				best_sol[bidx + pidx] = &subsol[pidx];

				auto& sol_info(best_sol[bidx + pidx]);
				sol_info->m_show_pop_id = 0;
				sol_info->m_pop_id = pidx;
				sol_info->m_sol_id = 0;
				sol_info->m_pop_ratio = double(sol_info->m_pop_id) / double(m_sub_pops.size());
				sol_info->m_sol_id_ratio = 0;
			}

			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
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