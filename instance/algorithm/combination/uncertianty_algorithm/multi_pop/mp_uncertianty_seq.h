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


#ifndef MP_UNCERTIANTY_SEQ_H
#define MP_UNCERTIANTY_SEQ_H




#include "pop_mp_main_uncertianty_seq.h"
#include "pop_mp_sub_uncertianty_seq.h"
#include "../../../../../core/algorithm/multi_population.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"

#include<list>

#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/multi_point_distribution/buffer_multi_point_distribution.h>
#include <custom/buffer/algorithm/combination/item_multi_opts.h>
#endif //  OFEC_DEMO


// for test
#include "../../../combination/sequence/SP/sp_operator.h"
#include "../../../template/framework/uncertianty/distance_calculator.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"





namespace ofec {

	struct IdManager {
		std::list<int> m_ids;
		int m_max_id = 0;
		void initialize() {
			m_ids.clear();
			m_max_id = 0;
		}
		int getId() {
			if (m_ids.empty()) {
				return m_max_id++;
			}
			else {
				int id = m_ids.back();
				m_ids.pop_back();
				return id;
			}
		}
		void insertId(int id) {
			m_ids.push_back(id);
		}
	};

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class MultiPopUncertiantySeq : public AlgorithmSeq {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename OpType::SolutionType;
			using InterpreterType = typename OpType::InterpreterType;
			using MainPopType = typename PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>;
			using SubPopType = typename PopulationMPSubUncertiantySeqEnsemble
				<TDistanceCalculator, TSequenceOperator>;

		protected:
			
			MainPopType m_main_pop;
			MultiPopulation<SubPopType> m_sub_pops;
			std::shared_ptr<EvaluationStrategyBase<SolutionType>> m_eval_strategy;
			IdManager m_ids;

			//std::vector<std::pair<int, int>, int> m_map_overlapping_info;


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
			virtual void appendPop(std::list<std::unique_ptr<SolutionType>>& inds);
			virtual void addNewlyPop();
			virtual void filterPop();
	};

	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopUncertiantySeq
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
		m_main_pop.resize(std::get<int>(
			m_param->at("population size")), 
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
		inline void MultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::addNewlyPop()
	{
		if (m_sub_pops.size() ==0) {
			m_main_pop.setSaveState(true);
			
			if (m_main_pop.judgeReady()) {
				std::list<std::unique_ptr<SolutionType>> inds;
				for (int idx(0); idx < m_main_pop.size(); ++idx) {
					if (m_main_pop.judgeNewlyPop(idx)) {
						inds.clear();
						m_main_pop.getPop(idx, inds);
						appendPop(inds);
					}
				}
				m_main_pop.setSaveState(false);
			}

		}

	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::filterPop()
	{


		bool erase_flag(false);
		int popx_id(0);
		auto& pops(m_sub_pops.getPops());
		int cur_id(0);
		int relation_ship(-1);
		std::vector<int> belong(pops.size(),-1);

		struct relationShipInfo {
			double m_x_fit;
			double m_y_fit;
			double m_thread_hold = 0;
			double m_dis = 0;
		};


		std::vector<std::vector<relationShipInfo>> belong_dis;
		//belong_dis.resize(pops.size(), -1);

		belong_dis.resize(pops.size());
		for (auto& it : belong_dis) {
			it.resize(pops.size());
		}

		//std::vector<double> fitness(pops.size(), -std::numeric_limits<double>::max());
		for (int idx(0); idx < pops.size(); ++idx) {
			for (int idy(0); idy < pops.size(); ++idy) {

				auto& info = belong_dis[idx][idy];

				info.m_x_fit = pops[idx]->getRepresentative().fitness();
				info.m_y_fit = pops[idy]->getRepresentative().fitness();
				info.m_thread_hold = pops[idx]->getDistanceCalculator()->radiusThreadhold();
				info.m_dis = pops[idx]->getDistanceCalculator()->disToPop(pops[idy]->getRepresentative());


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
			if (belong[idx]!=-1&&pops[belong[idx]]->id()!=-1) {
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
		inline void MultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::appendPop(std::list<std::unique_ptr<SolutionType>>& inds)
	{

		std::unique_ptr<SubPopType> curPop(new SubPopType);
		curPop->initParamenters(m_problem.get(), this);
		curPop->setEvalStrategy(m_eval_strategy);
		curPop->initialize(inds);
		curPop->setId(m_ids.getId());
		m_eval_strategy->clearMemory(curPop->id());
		curPop->updateMemory();
		m_sub_pops.append(curPop);
	}


	template <
		template<class> class TEvaluationStrategy,
		template<class> class TDistanceCalculator,
		class TAdaptor>
		inline void MultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::run_()
	{

		int iter(0);
		while (!terminating()) {
			++iter;
			//	std::cout << "alg id\t" << this << std::endl;
			m_main_pop.evolve(m_problem.get(), this, m_random.get());
			m_main_pop.updateRrepresentative();

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].evolve(m_problem.get(), this, m_random.get());
			}
			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateRrepresentative();
			}

#ifdef  OFEC_DEMO
			updateBuffer();
#endif

			addNewlyPop();
			filterPop();
			for (auto it(m_sub_pops.begin()); it != m_sub_pops.end(); ) {
				if ((**it).judgeConverged()) {
					m_main_pop.insertBestFitness((**it));
					it = m_sub_pops.remove(it);
				}
				else ++it;
			}


			m_main_pop.updateSurvivors();
			m_main_pop.updateMemory();
			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				m_sub_pops[idx].updateSurvivors();
				m_sub_pops[idx].updateMemory();
			}
#ifdef  OFEC_DEMO
			updateBuffer();

		//	m_problem->updateLocalOptimaInfo();
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
		inline void MultiPopUncertiantySeq
		<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::updateBuffer() {

		m_main_pop.updateIndisInfo();
		for (int idx(0); idx < m_sub_pops.size(); ++idx) {
			m_sub_pops[idx].updateIndisInfo();
		}


		ofec_demo::BufferMultiPointDistribution::m_numOpts
			= m_problem->optBase().numberVariables();
		ofec_demo::BufferMultiPointDistribution::m_numPops
			= m_sub_pops.size() + 1;


		int numPops(m_sub_pops.size());

		using sol_info_type = ofec_demo::BufferMultiPointDistribution::SolDataInfo;



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


			{
				int idx = 0;
				pop = &m_main_pop;
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
			using pop_type = PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>;
			std::vector<const pop_type* > sub_pops(m_sub_pops.size());

			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				sub_pops[idx] = &m_sub_pops[idx];
			}
			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
				drawMultiPop(sub_pops, m_problem.get());


		}


		{
			
			auto& pop_sols(ofec_demo::ItemMultiPoints::m_pop_sols);
			auto& sol_color(ofec_demo::ItemMultiPoints::m_pop_sol_fitcolor);
			pop_sols.clear();
			sol_color.clear();

			using pop_type = PopulationUncertiantySeqEnsemble
				<TDistanceCalculator, TAdaptor>;
			std::vector<const pop_type*> sub_pops(m_sub_pops.size()+1);
			sub_pops.front() = &m_main_pop;
			for (int idx(0); idx < m_sub_pops.size(); ++idx) {
				sub_pops[idx+1] = &m_sub_pops[idx];
			}

			for (auto& pop : sub_pops) {
				for (int solid(0); solid < pop->size(); ++solid) {
					auto& cur_sol((*pop)[solid]);
					const SolutionBase* a = &cur_sol;
					pop_sols.push_back(a);
					
				    double obj= cur_sol.objective()[0];
					obj =mapRealInside<double>
						(obj, obj_range.first, obj_range.second,
							0, 1);


					sol_color.push_back(obj);
				}
			}





		}

		//{
		//	using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;
		//	SolDataInfo cur;

		//	auto& solInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_sols;
		//	auto& obestInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_obest;
		//	auto& optInfo = ofec_demo::BufferMultiPointDistribution::m_subpop_opt;

		//	solInfo.clear();
		//	obestInfo.clear();
		//	optInfo.clear();

		//	const PopulationUncertiantySeqEnsemble
		//		<TDistanceCalculator, TAdaptor>* pop = nullptr;
		//	auto& opt_sol(m_his_opts.front());
		//	std::vector<int> fit_order;

		//	for (int idpop(0); idpop < m_sub_pops.size(); ++idpop) {
		//		cur.m_show_pop_id = idpop + 1;
		//		cur.m_opt_id = 0;
		//		cur.m_pop_id = idpop+1;
		//		cur.m_pop_ratio = double(cur.m_pop_id) / double(m_sub_pops.size());
		//		pop = &m_sub_pops[idpop];
		//		{
		//			std::vector<int> order;
		//			order.resize(pop->size());
		//			for (int idx(0); idx < pop->size(); ++idx) {
		//				order[idx] = idx;
		//			}
		//			std::sort(order.begin(), order.end(),
		//				[&](int a, int b) {
		//				return (*pop)[a].fitness() > (*pop)[b].fitness();
		//			});

		//			fit_order.resize(order.size());
		//			for (int idx(0); idx < fit_order.size(); ++idx) {
		//				fit_order[order[idx]] = idx;
		//			}
		//		}

		//		for (int solid(0); solid < pop->size(); ++solid) {
		//			cur.m_sol_id = solid;
		//			cur.m_sol_id_ratio = double(solid) / pop->size();
		//			cur.m_fitness_order_ratio = double(fit_order[solid]) / pop->size();
		//			auto& sol((*pop)[solid]);
		//			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//				updateDisInfo(cur, pop, sol, obj_range, opt_sol, m_problem.get());
		//			solInfo.push_back(cur);
		//		}

		//		for (int popidy(0); popidy < m_sub_pops.size(); ++popidy) {
		//			if (popidy != idpop) {
		//				cur.m_pop_id = popidy+1;
		//				cur.m_pop_ratio = double(cur.m_pop_id)/m_sub_pops.size();
		//				cur.m_sol_id = 0;
		//				cur.m_sol_id_ratio = 0;
		//				cur.m_fitness_order_ratio = 0;
		//				auto& sol(m_sub_pops[popidy].getRepresentative());
		//				PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//					updateDisInfo(cur, pop, sol, obj_range, opt_sol, m_problem.get());
		//				obestInfo.push_back(cur);
		//			}
		//		}
		//		cur.m_pop_id = -1;
		//		cur.m_fitness_order_ratio = 0;
		//		for (auto solIter(m_his_opts.begin()); solIter != m_his_opts.end(); ++solIter) {
		//			cur.m_sol_id = solIter->id();
		//			cur.m_sol_id_ratio = double(cur.m_sol_id) / double(m_his_opts.size());
		//			auto& sol(*solIter);

		//			PopulationUncertiantySeqEnsemble<TDistanceCalculator, TAdaptor>::
		//				updateDisInfo(cur, pop, sol, obj_range, opt_sol, m_problem.get());
		//			optInfo.push_back(cur);
		//		}
		//	}
		//}


		{
			using SolDataInfo = ofec_demo::BufferMultiPointDistribution::SolDataInfo;
			std::vector<SolDataInfo*> best_sol(m_main_pop.getSurvivors().size() + m_sub_pops.size());
			std::vector<const SolutionType*> best_indis(best_sol.size());

			auto& maintsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_distr);
			auto& subsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_subpop);
			auto& optsol(ofec_demo::BufferMultiPointDistribution::m_best_sol_opt);
			{
				int bidx(0);	
				maintsol.resize(m_main_pop.getSurvivors().size());
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


		{
		
			
		}


		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

	//template <
	//	template<class> class TEvaluationStrategy,
	//	template<class> class TDistanceCalculator,
	//	class TAdaptor>
	//	inline void GeneticLearningUncertiantySeq
	//	<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::updateOptPointInfo() {

	//	auto& alg(this);
	//	auto& pro(m_problem.get());


	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
	//	auto& optBase = pro->optBase();
	//	auto& optSol = pro->optBase().variable(0);
	//	auto& pointDensity = ofec_demo::BufferPointDistribution::m_opt_infos;
	//	pointDensity.resize(optBase.numberVariables());
	//	for (int idx(0); idx < optBase.numberVariables(); ++idx) {
	//		//interpreter.updateEdges()
	//		auto& it(pointDensity[idx]);
	//		SolutionType sol(optBase.variable(idx));
	//		sol.evaluate(pro, alg, false);
	//		interpreter.updateEdges(pro, sol);
	//		it.m_dis2opt = sol.variableDistance(optSol, pro);
	//		it.m_fitness = 1.0;
	//		it.m_id = idx;
	//		double fit = 1.0;
	//		it.m_dis2pop = m_pop.getDistanceCalculator()->disToPop(sol);
	//	}

	//}
	//template <
	//	template<class> class TEvaluationStrategy,
	//	template<class> class TDistanceCalculator,
	//	class TAdaptor>
	//	inline void GeneticLearningUncertiantySeq
	//	<TEvaluationStrategy, TDistanceCalculator, TAdaptor>::updatePointDensity() {

	//	auto& alg(this);
	//	auto& pro(m_problem.get());

	//	auto& pop(m_pop);

	//	using namespace ofec_demo;
	//	auto& optSol = pro->optBase().variable(0);
	//	double dis(optSol.variableDistance(optSol, pro));
	//	auto& pointDensity = ofec_demo::BufferPointDistribution::m_cur_sol_infos;
	//	pointDensity.resize(pop.size());

	//	auto& inter(GET_ASeq(alg).getInterpreter<InterpreterType>());
	//	for (int idx(0); idx < pop.size(); ++idx) inter.updateEdges(pro, pop[idx]);

	//	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
	//	//	auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>);
	//	auto obj_range(pro->optBase().objRange());
	//	for (int idx(0); idx < pop.size(); ++idx) {

	//		//interpreter.updateEdges()
	//		auto& it(pointDensity[idx]);
	//		it.m_dis2opt = pop[idx].variableDistance(optSol, pro);
	//		it.m_fitness = pop[idx].objective()[0];
	//		if (it.m_fitness < obj_range.first) it.m_fitness = obj_range.first;
	//		else if (it.m_fitness > obj_range.second) it.m_fitness = obj_range.second;
	//		it.m_fitness = mapReal<double>(it.m_fitness, obj_range.first, obj_range.second, 0, 1);
	//		if (pro->optimizeMode()[0] == OptimizeMode::kMinimize) {
	//			it.m_fitness = 1.0 - it.m_fitness;
	//		}
	//		it.m_id = idx;
	//		it.m_dis2pop = m_pop.getDistanceCalculator()->disToPop(pop[idx]);
	//	}
	//}
#endif
}


#endif