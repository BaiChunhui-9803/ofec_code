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


#ifndef AMP_UNCERTIANTY_SEQ_H
#define AMP_UNCERTIANTY_SEQ_H




#include "pop_mp_main_uncertianty_seq.h"
#include "pop_mp_sub_uncertianty_seq.h"
#include "mp_uncertianty_seq.h"


#include "../../../../../core/algorithm/multi_population.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"


#include "../../../../../utility/clustering/fdc.h"

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
		class AdaptiveMultiPopUncertiantySeq : public AlgorithmSeq {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename OpType::SolutionType;
			using InterpreterType = typename OpType::InterpreterType;
			using MainPopType = typename PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>;
			using SubPopType = typename PopulationMPSubUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>;

		protected:

			
			//std::vector<std::unique_ptr<SolutionType>> m_temp_indis;
			Population<SolutionType> m_main_pop;
			MultiPopulation<SubPopType> m_sub_pops;
			std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_strategy;
			IdManager m_ids;

			int m_new_popNum = 200;
			int m_init_num = 50;
			int m_min_popNum = 20;
			//double m_radius = 0.85;


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
		inline void AdaptiveMultiPopUncertiantySeq
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
	//	m_main_pop.initParamenters(m_problem.get(), this);
	//	m_main_pop.setEvalStrategy(m_eval_strategy);
	//	m_main_pop.setId(m_ids.getId());

	//	m_main_pop.initialize(m_problem.get(), m_random.get());
	//	m_main_pop.updateMemory();



		addNewlyPop();

		for (int idx(0); idx < m_sub_pops.size(); ++idx) {
			m_sub_pops[idx].updateSurvivors();
			m_sub_pops[idx].updateMemory();
		}
	}

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void AdaptiveMultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::addNewlyPop()
	{
		int origin_size(m_main_pop.size());
		int cur_size = m_new_popNum;
		m_main_pop.resize(cur_size, m_problem.get(), m_problem->numberVariables());

		auto& inds(m_main_pop);
		for (int i = origin_size; i < inds.size(); ++i) {
			inds[i]->initialize(i, m_problem.get(), m_random.get());
			m_eval_strategy->evaluate(*inds[i], -1);
		}
		//std::vector<int> fitness_order;
		//std::vector<int> belong_idx;
		//std::vector<int> distance_order;
		//sortByFitnessDis(inds, fitness_order, belong_idx, distance_order, m_problem.get());

		//std::vector<int> pop_num(belong_idx.size(), 0);
		//std::vector<double> min_dis(belong_idx.size(), 0);
		//
		//for (auto& it : fitness_order) {
		//	min_dis[it] = inds[it]->variableDistance(*inds[belong_idx[it]], m_problem.get());
		//	if (inds[it]->variableDistance(*inds[belong_idx[it]],m_problem.get()) < m_radius) {
		//		belong_idx[it] = belong_idx[belong_idx[it]];
		//	}
		//	else{
		//		belong_idx[it] = it;
		//	}
		//	++pop_num[belong_idx[it]];
		//}
		//
		//std::vector<int> pop_ids;
		//for (int idx(0); idx < pop_num.size(); ++idx) {
		//	if (pop_num[idx] >= m_init_num) {
		//		pop_ids.push_back(idx);
		//	}
		//}
		//
		FDC<SolutionType> clu(m_main_pop, m_problem.get());
		clu.calFitDist();
		///clu.calculate_objective_distance();
		clu.clusterBySize(m_min_popNum);
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

		
		int curIdx(0);
		for (int idx(0); idx < inds.size(); ++idx) {
			if (inds[idx] != nullptr&&curIdx!=idx) {
				inds[curIdx++] = std::move(inds[idx]);
			}
		}
	
		m_main_pop.resize(curIdx, m_problem.get(), m_problem->numberVariables());
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void AdaptiveMultiPopUncertiantySeq
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
								pops[belong[idy]]->getRepresentative().fitness()
								) {
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


		//std::vector<double> re_fit(m_sub_pops.size());
		//for (int idx(0); idx < m_sub_pops.size(); ++idx) {
		//	re_fit[idx] = m_sub_pops[idx].getRepresentative().fitness();
		//}
		//std::vector<int> fit_order(m_sub_pops.size());
		//for (int idx(0); idx < fit_order.size(); ++idx) {
		//	fit_order[idx] = idx;
		//}
		//std::sort(fit_order.begin(), fit_order.end(), [&](int a, int b) {
		//	return re_fit[a] < re_fit[b];
		//});

		//for (int re_idx(0); re_idx < fit_order.size(); ++re_idx) {
		//	int idx = fit_order[re_idx];
		//	if (pops[idx] != nullptr) {
		//		erase_flag = false;
		//		popx_id = pops[idx]->id();
		//		for (int re_idy(0); re_idy < re_idx; ++re_idy) {
		//			int idy = fit_order[re_idy];
		//			if (pops[idy] != nullptr) {
		//				int ret = pops[idx]->judgeRelationship(*pops[idy]);
		//				if (ret > 0) {
		//					pops[idx]->mergePopByFitnessDis()
		//				}
		//				else if(ret<0){
		//					erase_flag = true;
		//				}
		//			}
		//		}
		//	}
		//}



		//for (int idx(0); idx < m_sub_pops.size(); ++idx) {
		//	if (pops[idx] != nullptr) {
		//		erase_flag = false;
		//		popx_id = pops[idx]->id();
		//		for (int idy(0); idy < m_sub_pops.size(); ++idy) {
		//			if (pops[idy] != nullptr) {
		//				int ret = pops[idx]->judgeRelationship(*pops[idy]);
		//				if (ret > 0) {
		//					pops[idx]
		//				}
		//				else {

		//				}
		//				
		//				auto& it2(pops[idy]);
		//				//pops[idx]->mergePopByFitnessDis(*pops[idy]);
		//			//	pops[idx]->mergePopByFitnessDis(pops[idy]);
		//				if (pops[idy]->judgeSearchInside(*pops[idx])) {
		//					pops[idx]->mergePopByFitnessDis(*pops[idy]);
		//					pops[idy].reset(nullptr);
		//					erase_flag = true;
		//				}
		//				else if (m_sub_pops[idx].judgeSearchAreaInside(m_sub_pops[idy])) {
		//					pops[idx].reset(nullptr);
		//					erase_flag = true;
		//				}
		//			}

		//			if (erase_flag) {
		//				pops[idy].reset(nullptr);
		//				m_ids.insertId(pop_id);
		//			}
		//		}
		//	}
		//}
	}

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void AdaptiveMultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::
		appendPop(std::vector<std::unique_ptr<SolutionType>>& inds)
	{

		std::unique_ptr<SubPopType> curPop(new SubPopType);
		curPop->initParamenters(m_problem.get(), this);
		curPop->setEvalStrategy(m_eval_strategy);
		curPop->initialize(inds);
		curPop->setId(m_ids.getId());
		m_eval_strategy->clearMemory(curPop->id());
	//	curPop->updateMemory();
		m_sub_pops.append(curPop);
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void AdaptiveMultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::run_()
	{

		int iter(0);
		while (!terminating()) {
			++iter;
			//	std::cout << "alg id\t" << this << std::endl;
			//m_main_pop.evolve(m_problem.get(), this, m_random.get());
			//m_main_pop.updateMemory();


			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].evolve(m_problem.get(), this, m_random.get());
			}

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateRrepresentative();
			}
			//addNewlyPop();
			filterPop();

			for (auto it(m_sub_pops.begin()); it != m_sub_pops.end(); ) {
				if ((**it).judgeConverged()) {
					m_main_pop.append((*it)->getRepresentative());
					it = m_sub_pops.remove(it);
				}
				else ++it;
			}

			if (m_sub_pops.size() <=1) {
				addNewlyPop();
			}

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateSurvivors();
				m_sub_pops[idx].updateMemory();
			}


#ifdef  OFEC_DEMO
			updateBuffer();
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
		inline void AdaptiveMultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::updateBuffer() {


		for (int idx(0); idx < m_sub_pops.size(); ++idx) {
			m_sub_pops[idx].updateIndisInfo();
		}


		ofec_demo::BufferMultiPointDistribution::m_numOpts
			= m_problem->optBase().numberVariables();
		ofec_demo::BufferMultiPointDistribution::m_numPops
			= m_sub_pops.size();
		int numPops(m_sub_pops.size());

		using sol_info_type = ofec_demo::BufferMultiPointDistribution::SolDataInfo;

		auto& pop_idxs(ofec_demo::BufferMultiPointDistribution::m_popIdxs);
		pop_idxs.resize(numPops);
		for (int idx(0); idx < pop_idxs.size(); ++idx) {
			pop_idxs[idx] = m_sub_pops[idx].id();
		}

		auto& popDivLines(ofec_demo::BufferMultiPointDistribution::m_pop_div_val);
		popDivLines.resize(numPops);
		for (int idx(0); idx < pop_idxs.size(); ++idx) {
			popDivLines[idx] = m_sub_pops[idx].getDistanceCalculator()->normalize01(m_sub_pops[idx].getDistanceCalculator()->radiusThreadhold());
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



		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
			updateHisOpts(m_problem.get(), this);

		auto& obj_range(m_problem->optBase().objRange());
		{

			auto& sol_info(ofec_demo::BufferMultiPointDistribution::m_cur_sols);
			sol_info.clear();
			sol_info_type cur_sol_info;
			auto& opt_sol(m_his_opts.front());
			const PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>* pop = nullptr;
			for (int idx(0); idx < 1; ++idx) {
				pop = &m_sub_pops[idx];
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
		//		pop = &m_sub_pops[idx - 1];
		//		cur_sol_info.m_show_pop_id = idx;
		//		cur_sol_info.m_opt_id = 0;

		//		for (int idy(1); idy < pop_idxs.size(); ++idy) {
		//			cur_sol_info.m_pop_id = idy;
		//			if (idx != idy) {
		//				pop2 = &m_sub_pops[idy - 1];
		//				/*for (int solid(0); solid < pop2->size(); ++solid)*/ {
		//					cur_sol_info.m_sol_id = 0;
		//					//auto& sol((*pop2)[solid]);
		//					auto& sol((*pop2).getRepresentative());
		//					PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//						updateDisInfo(cur_sol_info, pop, sol, obj_range, opt_sol, m_problem.get());


		//					//cur_sol_info.m_dis2pop = pop->getDistanceCalculator()->disToPop(sol);
		//					//cur_sol_info.m_dis2opt = opt_sol.variableDistance(sol, m_problem.get());
		//					//cur_sol_info.m_objectives = mapReal<double>
		//					//	(sol.objective()[0], obj_range.first, obj_range.second,
		//					//		0, 1);
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
				<TDistanceCalculator, TAdaptor>* pop = nullptr;
			for(int pidx(0);pidx<m_sub_pops.size();++pidx)
			{
				pop = &m_sub_pops[pidx];
				cur_sol_info.m_show_pop_id = pidx;
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
			using pop_type = PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>;
			std::vector<const pop_type* > sub_pops(m_sub_pops.size());

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				sub_pops[idx] = &m_sub_pops[idx];
			}
			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
				drawMultiPop(sub_pops, m_problem.get(),0);


		}


		//{
		//	std::vector<const SolutionType*> best_indis(m_main_pop.size() + m_sub_pops.size());

		//	auto& best_sol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
		//	best_sol.resize(best_indis.size());
		//	{


		//		int bidx(0);
		//		for (; bidx < m_main_pop.size(); ++bidx) {
		//			best_indis[bidx] = &m_main_pop[bidx];
		//			best_sol[bidx].m_show_pop_id = -1;
		//			best_sol[bidx].m_pop_id = 0;
		//			best_sol[bidx].m_opt_id = 0;

		//		}
		//		for (int pidx(0); pidx < m_sub_pops.size(); ++pidx) {
		//			best_indis[bidx + pidx] = &m_sub_pops[pidx].getRepresentative();
		//			best_sol[bidx + pidx].m_show_pop_id = -1;
		//			best_sol[bidx + pidx].m_pop_id = pidx + 1;
		//			best_sol[bidx + pidx].m_opt_id = 0;
		//		}
		//		PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//			updateBestSolInfo(best_sol, best_indis, obj_range, m_problem.get());;

		//	}
		//}


		{
			using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;
			std::vector<SolDataInfo*> best_sol(m_main_pop.size() + m_sub_pops.size());
			std::vector<const SolutionType*> best_indis(best_sol.size());

			auto& maintsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
			auto& subsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_subpop);
			auto& optsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_opt);
			{
				int bidx(0);
				maintsol.resize(m_main_pop.size());
				for (; bidx < m_main_pop.size(); ++bidx) {
					best_indis[bidx] = &m_main_pop[bidx];
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