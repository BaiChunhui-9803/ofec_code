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

#ifndef OFEC_POP_MP_MAIN_UNCERTIANTY_SEQ_ENSEMBLE_H
#define OFEC_POP_MP_MAIN_UNCERTIANTY_SEQ_ENSEMBLE_H


#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

#include "../pop_uncertainty_seq_ensemble.h"
#include <algorithm>
#include <memory>


// for test
#include "../../../combination/sequence/SP/sp_operator.h"
#include "../../../template/framework/uncertianty/distance_calculator.h"
#include "../../../template/framework/uncertianty/evaluation_strategy.h"



namespace ofec {
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		class PopulationMPMainUncertiantySeqEnsemble :
		virtual public PopulationUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator> {
		public:
			using OpType = typename TSequenceOperator;
			using SolutionType = typename TSequenceOperator::SolutionType;
			using InterpreterType = typename TSequenceOperator::InterpreterType;

		public:
			virtual void setSolution(int id, const SolutionType& indi)override
			{
				int his_id = m_individuals[id]->id();
				PopulationUncertiantySeqEnsemble::setSolution(id, indi);				
				m_individuals[id]->setId(his_id);
			}
			virtual void insert(std::vector<std::unique_ptr<SolutionType>>& pops);
		
			template<typename ... Args>
			void resize(size_t size, Problem *pro, Args&& ...args){
				int origin_size(m_individuals.size());
				if (origin_size > size) {
					for (int idx(size); idx < origin_size; ++idx) {
						pushIds(m_individuals[idx]->id());
					}
				}
				PopulationUncertiantySeqEnsemble::resize(size, pro, std::forward<Args>(args)...);
				if (size > origin_size) {
					for (int idx(origin_size); idx < size; ++idx) {
						m_individuals[idx]->setId(assignIds());
					}
				}
			}
			virtual void deleteIndisFit(int numPop) {
				std::sort(m_individuals.begin(), m_individuals.end(), [](
					const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b
					) {
					return a->fitness() > b->fitness();
				});
				resize(numPop, m_problem.get());
			}


			//virtual int selectSurvivors()override {
			//	int origin_size(m_individuals.size());
			//	resize(m_individuals.size() + m_offsprings.size() , m_problem.get(), m_problem->numberVariables());
			//	for (int idx(origin_size); idx < m_individuals.size(); ++idx) {
			//		setSolution(idx, m_offsprings[idx - origin_size]);
			//	}
			//	reduceByFitness2(origin_size);
			//	if (m_random->uniform.next() > m_global_learn_pro) {
			//		m_type = 1;
			//	}
			//	else m_type = 0;
			//	return 0;
			//}

			virtual void insertBestFitness(Population<SolutionType> &pop) {
				int origin_size(m_individuals.size());
				resize(m_individuals.size() + pop.size(), m_problem.get(), m_problem->numberVariables());

				for (int idx(origin_size); idx < m_individuals.size(); ++idx) {
					setSolution(idx, pop[idx - origin_size]);
				}
				reduceByFitness2(origin_size);
			}
			virtual void reduceByFitness2(int maxPop) {
				auto& pops(m_individuals);
				if (pops.size() <= m_init_popsize)return;
				std::vector<int> fitness_order;
				std::vector<int> belong_idx;
				std::vector<int> dis_order;
				sortByFitnessDis(pops, fitness_order, belong_idx, dis_order, m_problem.get());
				for (int idx(pops.size() - 1); idx >= maxPop; --idx) {
					int cur_id(dis_order[idx]);
					pushIds(m_individuals[cur_id]->id());
					m_individuals[cur_id].reset(nullptr);
					
				}
				int cur_size(0);
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					if (m_individuals[idx] != nullptr) {
						if (idx != cur_size) {
							m_individuals[cur_size] = std::move(m_individuals[idx]);
							++cur_size;
						}
					}
				}
				PopulationUncertiantySeqEnsemble::resize(cur_size, m_problem.get(), m_problem->numberVariables());

			}


			virtual int selectSurvivors() override {
				if (m_type == 0) {
					auto& pop(m_individuals);
					m_num_improve = 0;
					for (int idx(0); idx < pop.size(); ++idx) {
						if (m_offsprings[idx].fitness() >= pop[idx]->fitness()) {
							
							setSolution(idx, m_offsprings[idx]);
							pop[idx]->setImproved(true);
							++m_num_improve;
							m_stagnated_times[idx] = 0;
						}
						else {
							pop[idx]->setImproved(false);
							m_his[pop[idx]->id()]->emplace_back( new SolutionType(m_offsprings[idx]));
							m_stagnated_times[pop[idx]->id()]++;
						}
					}
					mergeByAttraction();
				}
				else {
					int origin_size(m_individuals.size());
					resize(m_individuals.size() + m_offsprings.size(), m_problem.get(), m_problem->numberVariables());

					for (int idx(origin_size); idx < m_individuals.size(); ++idx) {
						setSolution(idx, m_offsprings[idx - origin_size]);
					}
					reduceByFitness2(origin_size);
					//std::vector<SolutionType> newly_inds(std::move(m_offsprings));
					//int origin_size(m_individuals.size());
					//resize(m_individuals.size() + newly_inds.size(), m_problem.get(),m_problem->numberVariables());
					//for (int idx(0); idx < newly_inds.size(); ++idx) {
					//	setSolution(idx + origin_size, newly_inds[idx]);
					//	m_individuals[idx]->setImproved(true);
					//	m_stagnated_times[idx] = 0;
					//}
					//for (int idx(0); idx < origin_size; ++idx) {
					//	m_individuals[idx]->setImproved(false);
					//	m_stagnated_times[idx]++;
					//}
				}

				if (m_save_state) {
					m_type = 0;
				}
				else {
					if (m_random->uniform.next() > m_global_learn_pro) {
						m_type = 1;
					}
					else m_type = 0;
				}


				//if (m_individuals.size() < m_init_size) {
				//	m_type = 1;
				//	resizeOffsprings(m_init_size - m_individuals.size());
				//}
				//else {
				//	if (m_type == 1) {
				//		m_type = 0;
				//		resizeOffsprings(m_individuals.size());
				//	}
				//}
				//modifyStagnatedTimes();
				return 0;
				
			}

			virtual void setEvolvoType(int type) {
				m_type = type;
			}




			virtual void initialize(Problem *pro, Random *rnd)override {
				PopulationUncertiantySeqEnsemble::initialize(pro, rnd);
				m_type = 0;
				for (int idx(0); idx < m_individuals.size(); ++idx) {
					m_individuals[idx]->setId(assignIds());
				}
			}
			virtual bool judgeNewlyPop(int idIndi) {
				return m_individuals[idIndi]->isImproved()
					&& m_his[m_individuals[idIndi]->id()]->size() >= m_pop_init_size;
			}

			virtual void getPop(int idIndi,std::list<std::unique_ptr<SolutionType>>& inds) {
				inds.clear();
				auto & cur_list(*m_his[m_individuals[idIndi]->id()]);
				inds.splice(inds.begin(), cur_list);
				inds.emplace_back(new SolutionType(*m_individuals[idIndi]));
			}

			virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override {
				//return 0;
				m_actions.setPro(m_action_pro[m_type]);
				return PopulationUncertiantySeqEnsemble::evolve(pro, alg, rnd);
			}

			virtual bool judgeReady() {
				for (auto& it : m_individuals) {
					if (m_his[it->id()]->size() < m_pop_init_size) {
						return false;
					}
				}
				return true;
			}

			void setSaveState(bool flag) {
				m_save_state = flag;
			}

	    protected:
			virtual bool judgeInside(SolutionType& cur, SolutionType& other);
			virtual void mergeInds(int cur_id, int other);



			virtual void pushIds(int id) {
				m_his[id]->clear();
				m_stagnated_times[id] = m_init_stagnated_times;
				m_his_ids.push_back(id);
			}

			virtual int assignIds() {
				if (!m_his_ids.empty()) {
					int id(m_his_ids.back());
					m_his_ids.pop_back();
					m_stagnated_times[id]= m_init_stagnated_times;
					m_his[id]->clear();
					return id;
				}
				else {
					int id = m_his.size();
					m_his.emplace_back(new std::list<std::unique_ptr<SolutionType>>());
					m_stagnated_times.push_back(m_init_stagnated_times);
					return id;
				}
			}
			////virtual void resizePop();
			virtual void mergeByAttraction();
			virtual void modifyStagnatedTimes() {
				for (auto& it : m_stagnated_times) it = std::min(it, m_max_stagnated_times);
			}



		protected:

			int m_type = 0;
			std::vector<std::vector<double>> m_action_pro =
			{ {0.2,0.8,0} ,{1.0,0,0} };

			double m_global_learn_pro = 0.8;


			bool m_save_state = false;

			//std::vector<std::unique_ptr<Solution>> m_his;

			//std::vector<double> m_action_pro_0 = { 0.2,0.8,0 };
			//std::vector<double> m_action_pro_1 = {1.0,0,0};

			int m_init_stagnated_times = 2000;
			int m_max_stagnated_times = 2e5 + 5;
			int m_max_size = 80;
			int m_init_size = 50;
			int m_pop_init_size = 50;
			double m_inner_radius = 0.2;


			std::vector<int> m_stagnated_times;
			std::vector<std::unique_ptr<std::list<std::unique_ptr<SolutionType>>>> m_his;
			std::list<int> m_his_ids;
	};

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline void PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::
		insert(std::vector<std::unique_ptr<SolutionType>>& pops)
	{
		resize(m_individuals.size()+1, m_problem.get(), m_problem->numberVariables());
		int best_id(0);
		for (int idx(0); idx < pops.size(); ++idx) {
			if (pops[best_id]->fitness() < pops[idx]->fitness()) {
				best_id = idx;
			}
		}
		auto& list(m_his[m_individuals.back()->id()]);
		for (int idx(0); idx < pops.size(); ++idx) {
			if (idx != best_id) {
				list->emplace_back(std::move(pops[idx]));
			}
		}
		setSolution(m_individuals.size() - 1, *pops[best_id]);
		//mergeByAttraction();
	}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
		inline bool PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::judgeInside(SolutionType& cur, SolutionType& other)
	{
		if (cur.fitness() > other.fitness() && cur.variableDistance(other,m_problem.get()) <= m_inner_radius) {
			return true;
		}
		else {
			for (auto& it : *m_his[cur.id()]) {
				if (it->fitness() > other.fitness() && it->variableDistance(other,m_problem.get())<=m_inner_radius) {
					return true;
				}
			}
			return false;
		}
	}
	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	inline void PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::mergeInds(int cur_id, int other_id)
	{
		auto& cur(*m_individuals[cur_id]);
		auto& other(*m_individuals[other_id]);
		m_his[cur.id()]->splice(m_his[cur.id()]->begin(),*m_his[other.id()]);
		pushIds(other.id());
		m_his[cur.id()]->emplace_back(m_individuals[other_id].release());
		m_individuals[other_id].reset(nullptr);
	}



	//template <
	//	template<class> class TDistanceCalculator,
	//	class TSequenceOperator>
	//	inline void PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>::resizePop()
	//{

	//}

	template <
		template<class> class TDistanceCalculator,
		class TSequenceOperator>
	inline void PopulationMPMainUncertiantySeqEnsemble<TDistanceCalculator, TSequenceOperator>
		::mergeByAttraction()
	{
		auto& pops(m_individuals);
		std::vector<int> fitness_order(m_individuals.size(), -1);
		for (int idx(0); idx < fitness_order.size(); ++idx) {
			fitness_order[idx] = idx;
		}
		std::sort(fitness_order.begin(), fitness_order.end(),
			[&](int a,int b) {
			return m_individuals[a]->fitness() < m_individuals[b]->fitness();
		});

		for (int idx(0); idx < fitness_order.size(); ++idx) {
			int other_id(fitness_order[idx]);
			for (int idy(idx + 1); idy < fitness_order.size(); ++idy) {
				int cur_id(fitness_order[idy]);
				if (judgeInside(*m_individuals[cur_id], *m_individuals[other_id])) {
					{
						double fit1(m_individuals[cur_id]->fitness());
						double fit2(m_individuals[other_id]->fitness());
						double dis(m_individuals[cur_id]->variableDistance(*m_individuals[other_id], m_problem.get()));
					}
					mergeInds(cur_id, other_id);
					break;
				}
			}
		}

		//std::vector<int> belong_idx;
		//std::vector<int> dis_order;
		//sortByFitnessDis(pops, fitness_order, belong_idx, dis_order,m_problem.get());
		//
		//for (auto idIter(dis_order.rbegin()); idIter != dis_order.rend(); ++idIter) {
		//	int other_id(*idIter);
		//	int cur_id(belong_idx[other_id]);

		//	if (judgeInside(*m_individuals[cur_id], *m_individuals[other_id])) {

		//		{
		//			double fit1(m_individuals[cur_id]->fitness());
		//			double fit2(m_individuals[other_id]->fitness());
		//			double dis(m_individuals[cur_id]->variableDistance(*m_individuals[other_id], m_problem.get()));
		//		}
		//		mergeInds(cur_id, other_id);
		//	}
		//}
		//
		int cur_size(0);
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			if (m_individuals[idx] != nullptr ) {
				if (idx != cur_size)
					m_individuals[cur_size] = std::move(m_individuals[idx]);
				++cur_size;
			}
		}
		m_individuals.resize(cur_size);
		resize(cur_size, m_problem.get(), m_problem->numberVariables());
	}

}
#endif

