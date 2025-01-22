/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Jun 10, 2019 by Li Zhou

#ifndef OFEC_AMP_CC_CSIWDN_H
#define OFEC_AMP_CC_CSIWDN_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../dyn_sa_multipop.h"
#include "../cc_framework.h"
#include "../../../../../core/global.h"
#include <iomanip>      // std::setprecision

#define USE_LSTM

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../datum/datum_inclusion.h"




namespace ofec {
	template <typename TPop1, typename TPop2, typename TInd>
	class AMP_CC_CSIWDN : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(AMP_CC_CSIWDN)

		using PopulationType = CC_framework<TPop1, TPop2, TInd>;
		std::mutex mutex_AMP_out;

	public:

		void initialize_(Environment* env) override;
		void run_(Environment* env) override;
		void initMultiPop(Environment* env);
		const std::pair<int,int> & sourceStartEndIdx() { return m_source_idx; }

	protected:
		void addInputParameters();
		
		int evolve(Environment* env) ;
#ifdef OFEC_DEMO
		void updateBuffer();
#endif

		int competition(Environment* env);
		void computeNodeDistance(std::vector<Real> &vect, Environment* env);
		void addNewPop(size_t num, CSIWDN::InitType type, Environment* env);
		void output(Environment* env);
		bool isStableByPopNum();
		bool isStableByObj();
		void updateDistance(size_t idx, Environment* env);
		bool isStable();
		void setUseEpanet(Environment* env);
		std::vector<Real> calCoverageRate(Environment* env);
		void recordIntermediateResult(Environment* env);
		Real calStandardDeviation(const VarCSIWDN &var, const VarCSIWDN &opt_var, Environment* env);
		void record();
		void updateStartEndSouce(Environment* env);
	

		void initializePop(Environment* env);

		

	protected:

		DynSaMultiPop<PopulationType> m_subpop;

		int m_subpop_size;
		int m_num_add;
		std::string m_path_result;
		size_t m_data_update_eval;
		int m_iteration = 0;
		bool m_is_first_add = true;
		std::vector<int> m_num_be_removed;
		size_t m_num_curr_pop;
		size_t m_initial_num_pop;

		size_t m_k;
		CSIWDN::ClusterType m_cluster_type;
		size_t m_standard_stable;
		size_t m_reserve_stable;
		Real m_max_coverge_rate;

		std::vector<std::vector<Real>> m_coverage_result;
		std::vector<Solution<VarCSIWDN>> m_his_best;
		std::vector<int> m_curr_sourceidx;
		std::vector<size_t> m_num_pop_result;
		std::vector<Real> m_curr_best_obj;
		std::vector<size_t> m_curr_evals;
		std::vector<std::vector<size_t>> m_id_his_each_source;

	//	size_t m_maximum_evaluations;
		size_t m_initial_phase;
		size_t m_final_phase;
		int m_eval_time;  //延迟的优化时间
		std::pair<int, int> m_source_idx;

		size_t m_increment;
		size_t m_cur_contaminant;

		size_t m_lstm_phase;
		bool m_is_epanet_time = false;
		bool m_use_LSTM;
		std::vector<Solution<VarCSIWDN>> m_last_sol;

		size_t m_num_pop_last_iter = 0;
		int m_z = 0;
		
		//= m_num_curr_pop;
	};



	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::addInputParameters() {

		m_input_parameters.add("number of sub_populations", new RangedSizeT(m_initial_num_pop, 1, 500, 5));
		m_input_parameters.add("subpopulation size", new RangedInt(m_subpop_size, 1, 500, 10));
		//m_input_parameters.add("maximum evaluations", new RangedSizeT(m_maximum_evaluations, 10, 10000000, 200000));

		m_input_parameters.add("alpha", new RangedSizeT(m_increment, 0, 200, 2));
		m_input_parameters.add("beta", new RangedSizeT(m_standard_stable, 0, 200, 4));
		m_input_parameters.add("initial clustering", new Enumeration(m_cluster_type,
			{ "AHC","KMeans","NOCluster"}, CSIWDN::ClusterType::kAHC));

		m_input_parameters.add("time window", new RangedInt(m_eval_time, 1, 10000000, 86400));


		m_lstm_phase = 1;
		
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::initialize_(Environment* env) {
		Algorithm::initialize_(env);

		auto pro = env->problem();
		

	//	auto& v = *m_param;
		//m_lstm_phase = v.get<int>("use lstm proportion");
		//m_eval_time = v.get<int>("time window");
		//m_subpop_size = v.get<int>("subpopulation size");
		m_num_add = m_initial_num_pop / 2;
		m_k = m_initial_num_pop;
		m_num_curr_pop = m_initial_num_pop;
		//m_initial_num_pop = v.get<int>("number of subpopulations");
	//	m_cluster_type = v.has("initial clustering") ? static_cast<CSIWDN::ClusterType>(v.get<int>("initial clustering")) : CSIWDN::ClusterType::kAHC;
	//	m_standard_stable = v.get<int>("beta");
		m_reserve_stable = m_standard_stable;
	//	m_increment = v.get<int>("alpha");
//		m_maximum_evaluations = v.get<int>("maximum evaluations");
		m_initial_phase = CAST_CSIWDN(pro)->phase();
		m_final_phase = CAST_CSIWDN(pro)->totalPhase();
		m_data_update_eval = m_maximum_evaluations / (m_final_phase - m_initial_phase + 1);
		m_subpop.clear();
		CAST_CSIWDN(pro)->setClusterType(m_cluster_type);

		m_cur_contaminant = 0;
		m_source_idx.first = 0;
		m_source_idx.second = 0;

		m_use_LSTM = false;


	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::initMultiPop(Environment* env) {
		auto pro = env->problem();
		CAST_CSIWDN(pro)->setInitType(CSIWDN::InitType::kRandom);
		CAST_CSIWDN(pro)->clustering(m_k, m_random.get());
		for (size_t i = 0; i < m_num_curr_pop; ++i) {  // initialize subpopulation 
			auto subpop = std::make_unique<PopulationType>(m_subpop_size, env);
			CAST_CSIWDN(pro)->setPopIdentifier(i);
			subpop->m_first_pop->initialize(env, m_random.get());
			subpop->m_second_pop->initialize(env, m_random.get());
			subpop->fillEach(m_random.get(), env, PopulationType::FillType::random, m_source_idx);
			subpop->updateBest(env);
			m_subpop.append(subpop);
		}



	}	
	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::initializePop(Environment* env) {
		auto pro = env->problem();
		m_id_his_each_source.resize(CAST_CSIWDN(pro)->numberSource());
		initMultiPop(env);

#ifdef OFEC_DATUM_MULTI_POP_H

		g_multi_pop.bindMultiPopulations(m_subpop);
		env->algorithm()->datumUpdated(env, g_multi_pop);

#endif // OFEC_DATUM_MULTI_POP_H

	
		m_num_pop_last_iter = m_num_curr_pop;
		m_z = 0;
		CAST_CSIWDN(pro)->setAlgorithmStart();
	
	}
	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::run_(Environment* env) {

		initializePop(env);
		auto pro = env->problem();
		int tag = kNormalEval;



		while (!terminating() && tag != kTerminate) {
		//	int tag = kNormalEval;
			recordIntermediateResult(env);
			tag = evolve(env);

#ifdef OFEC_DATUM_MULTI_POP_H

			g_multi_pop.bindMultiPopulations(m_subpop);
			env->algorithm()->datumUpdated(env, g_multi_pop);

#endif // OFEC_DATUM_MULTI_POP_H
		}
		if (m_use_LSTM) {
			CAST_CSIWDN(pro)->useLstmModel(false);
		}
		for (size_t k = 0; k < m_subpop.size(); ++k) {
			if (m_use_LSTM) {
				m_subpop[k].evaluate(env);
				m_subpop[k].updateBest(env);
			}
			m_last_sol.push_back(dynamic_cast<Solution<VarCSIWDN>&>(*m_subpop[k].m_best));
		}
		//for (auto iter = m_candidates.begin(); iter != m_candidates.end(); iter++) {
		//	for (int i = 0; i < m_last_sol.size(); i++) {
		//		if (dynamic_cast<Solution<VarCSIWDN>&> (**iter).objective()[0] == m_last_sol[i].objective()[0])
		//			break;
		//		if(dynamic_cast<Solution<VarCSIWDN>&> (**iter).objective()[0] != m_last_sol[i].objective()[0] && i == m_last_sol.size()-1)
		//			m_last_sol.push_back(dynamic_cast<Solution<VarCSIWDN>&> (**iter));
		//	}
		//}
		output(env);
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd> ::setUseEpanet(Environment *env) {
		auto pro = env->problem();
		if (m_is_epanet_time) {
			CAST_CSIWDN(pro)->useLstmModel(false);
		}
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::recordIntermediateResult(Environment* env) {
		auto pro = env->problem();
		m_coverage_result.push_back(calCoverageRate(env));
		m_num_pop_result.push_back(m_subpop.size());

		Real best = m_subpop[0].m_best->objective()[0];
		size_t idx = 0;
		for (size_t i = 1; i < m_subpop.size(); ++i) {
			if (best > m_subpop[i].m_best->objective()[0]) {
				best = m_subpop[i].m_best->objective()[0];
				idx = i;
			}
		}

		for (size_t i = 0; i < m_subpop.size(); ++i) {
			for (size_t j = 0; j < CAST_CSIWDN(pro)->numberSource(); ++j) {
				std::cout << m_subpop[i].m_best->variable().location(j) << "  ";
			}
			std::cout << std::endl;
		}

		m_his_best.push_back(*m_subpop[idx].m_best);
		m_curr_best_obj.push_back(best);
		std::cout << "current best obj:" << best << "  " <<"num of pops:" << m_subpop.size() <<std::endl;
		m_curr_evals.push_back(m_evaluations);
		m_curr_sourceidx.push_back(m_source_idx.second);
		m_id_his_each_source[m_curr_sourceidx.back()].push_back(m_his_best.size() - 1);
		//std::cout << "current best obj:" << best << std::endl;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	std::vector<Real> AMP_CC_CSIWDN<TPop1, TPop2, TInd>::calCoverageRate(Environment *env) {
		auto pro = env->problem();
		int number_source = CAST_CSIWDN(pro)->numberSource();
		std::vector <std::vector<bool>> coverage(number_source, std::vector<bool>(CAST_CSIWDN(pro)->numberNode(), false));
		std::vector<int> count(number_source, 0);
		int z = m_source_idx.second, q = m_source_idx.first;
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			if (m_subpop[i].m_best) {
				for (size_t j = 0; j < m_subpop[i].m_first_pop->size(); ++j) {
					for (size_t k = q; k <= z; k++) {
						int index = (*(m_subpop[i].m_first_pop))[j].variable().index(k) - 1;
						if (!coverage[k][index]) { coverage[k][index] = true; ++count[k]; }
					}
				}
			}
		}
		std::vector<Real> count_coverage(number_source, 0);
		for (size_t i = q; i <= z; i++) {
			count_coverage[i] = (Real)count[i] / (Real)CAST_CSIWDN(pro)->numberNode();
		}
		return count_coverage;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	Real AMP_CC_CSIWDN<TPop1, TPop2, TInd>::calStandardDeviation(const VarCSIWDN& var, const VarCSIWDN& opt_var, Environment* env) {
		auto pro = env->problem();
		Real sum = 0;
		size_t total_size = 0;
		for (size_t i = 0; i < CAST_CSIWDN(pro)->numberSource(); i++) {
			if (var.multiplier(i).size() != opt_var.multiplier(i).size())
				return -1;
			size_t size = var.multiplier(i).size();
			for (size_t j = 0; j < size; j++)
				sum += pow(var.multiplier(i)[j] - opt_var.multiplier(i)[j], 2);
			total_size += size;
		}
		return sqrt(sum / (Real)total_size);
	}

	template<typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::record() {

		//m_coverage_result.push_back(calCoverageRate());
		//m_num_pop_result.push_back(m_subpop.size());
		//
		//std::cout << "eval:" << m_evaluations << "\t";
		//for (size_t i = 0; i < m_subpop.size(); ++i) {
		//	if (m_subpop[i].m_best) {
		//		std::cout << "subpop" << i << " best obj:";
		//		std::cout << m_subpop[i].m_best->objective()[0] << "   ";
		//		std::cout << "\t";
		//	}
		//}
		//std::cout << std::endl;

		//Real best = 11;
		//for (size_t i = 0; i < m_subpop.size(); ++i) {
		//	if (m_subpop[i].m_best) {
		//		if (best > m_subpop[i].m_best->objective()[0])
		//			best = m_subpop[i].m_best->objective()[0];
		//	}
		//}
		//m_curr_best_obj.push_back(best);
		//m_curr_evals.push_back(m_evaluations);

		/*auto &best_sol = dynamic_cast<Solution<VarCSIWDN>&>(*m_candidates.front());
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		entry.push_back(best_sol.objectiveDistance(CAST_CSIWDN(pro)->optima()->objective(0)));
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
		std::cout << "FEs: " << entry[0]  << "\tBest obj: " <<  std::setprecision(10) << entry[1] << "\tTime eval: " << best_sol.timeEvaluate() << std::endl;
		for (size_t i = 0; i < CAST_CSIWDN(pro)->numberSource(); ++i) {
			std::cout << "locations:" << best_sol.variable().location(i) << "\tmutilplier:" << std::setprecision(3);
			for(size_t j = 0; j < best_sol.variable().multiplier(i).size(); ++j){
				std::cout << best_sol.variable().multiplier(i)[j]<<"  ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;*/
	}

	template<typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::updateStartEndSouce(Environment* env)
	{
		auto pro = env->problem();
		int q = 0;
		long histry_time = CAST_CSIWDN(pro)->phase() * CAST_CSIWDN(pro)->intervalTimeStep() - m_eval_time;
		while (histry_time > 0 && (q + 1) < CAST_CSIWDN(pro)->numberSource() && histry_time >= (*(m_subpop[0].m_first_pop))[0].variable().startTime(q + 1)) {
			q++;
		}
		m_source_idx.first = q;
		int z = 0;
		while ((z + 1) < CAST_CSIWDN(pro)->numSource() && CAST_CSIWDN(pro)->phase() >= ((*(m_subpop[0].m_first_pop))[0].variable().startTime(z + 1) / CAST_CSIWDN(pro)->intervalTimeStep())) {
			z++;
		}
		m_source_idx.second = z;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	int AMP_CC_CSIWDN<TPop1, TPop2, TInd>::evolve(Environment *env) {

		auto pro = env->problem();


		const Real eta = 0.4;
		Real sum = 0.0;
		size_t num_pop = m_subpop.size();
		std::vector<bool> flag(num_pop, true);
		std::vector<Real> pops_obj(num_pop, 0.0);
		for (size_t i = 0; i < num_pop; ++i) {
			Real obj = m_subpop[i].m_best->objective()[0];
			//if (m_curr_best_obj.back() < 1) {
			//	if (obj > 1)
			//		obj = 1;
				if (obj < 0.0001)
					obj = 0.0001;
			//}
			pops_obj[i] = 1.0 / obj;
			sum += obj;
		}
		std::vector<Real> prob(num_pop);
		Real sum_pro = 0;
		for (size_t i = 0; i < num_pop; i++) {
			prob[i] = pops_obj[i];
			sum_pro += prob[i];
		}
		for (size_t i = 0; i < num_pop; i++)
			prob[i] = eta * num_pop * prob[i] / sum_pro;

		for (size_t i = 0; i < m_subpop.size(); ++i)
			if (m_random->uniform.next() < prob[i])
				m_subpop[i].evolve(env, m_random.get(), m_source_idx);
		++m_iteration;



		int num_remove = competition(env);
		m_num_be_removed.push_back(num_remove);
		m_num_curr_pop -= num_remove;
		if (isStable()) {
			std::cout << "Stagnation, new pops added." << std::endl;
			m_num_be_removed.clear();
			if (m_is_first_add) {
				addNewPop(m_num_add, CSIWDN::InitType::kBeVisited, env);
				m_is_first_add = false;
			}
			else {
				int differ = int(m_num_pop_last_iter) - int(m_num_curr_pop);
				if (differ < m_num_add)
					++m_num_add;
				else if (differ > m_num_add)
					--m_num_add;
				if (m_num_add == 0)
					m_num_add = 1;
				addNewPop(m_num_add, CSIWDN::InitType::kBeVisited, env);
			}
			m_num_pop_last_iter = m_num_curr_pop;
		}

		if (m_use_LSTM) {
			m_is_epanet_time = false;
			if (!CAST_CSIWDN(pro)->isUseLstmModel()) {
				CAST_CSIWDN(pro)->useLstmModel(true);
				for (size_t k = 0; k < m_subpop.size(); ++k) {
					m_subpop[k].evaluate(env);
					m_subpop[k].updateBest(env);
				}
			}
		}

		if (m_evaluations >= (CAST_CSIWDN(pro)->phase() - m_initial_phase + 1) * m_data_update_eval) {
			int current_phase = int(m_evaluations) / int(m_data_update_eval) + m_initial_phase;
			if (current_phase <= CAST_CSIWDN(pro)->totalPhase()) {
				CAST_CSIWDN(pro)->updatePhase(current_phase);
				updateStartEndSouce(env);
				if (m_use_LSTM) {
					if (current_phase % m_lstm_phase == 0)
						m_is_epanet_time = true;
				}
				//std::cout << "The current phase:" << current_phase << std::endl;
			}
			if (m_use_LSTM)
				setUseEpanet(env);
			//m_candidates.clear();
			for (size_t k = 0; k < m_subpop.size(); ++k) {
				m_subpop[k].evaluate(env);
				m_subpop[k].updateBest(env);
			}

			if ((m_z + 1) < CAST_CSIWDN(pro)->numSource() && CAST_CSIWDN(pro)->phase() >= ((*(m_subpop[0].m_first_pop))[0].variable().startTime(m_z + 1) / CAST_CSIWDN(pro)->intervalTimeStep()) && CAST_CSIWDN(pro)->phase() < CAST_CSIWDN(pro)->totalPhase()) {
				m_standard_stable = m_reserve_stable;
				addNewPop(m_initial_num_pop, CSIWDN::InitType::kCluster, env);
				m_z++;
				//std::cout << "The next source" << std::endl;
			}
			if (CAST_CSIWDN(pro)->phase() % m_increment == 0)
				m_standard_stable += 1;
		}

#ifdef OFEC_DEMO
		updateBuffer();
#endif // OFEC_DEMO


		return kNormalEval;
	}

#ifdef OFEC_DEMO
	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_subpop.size());
		for (size_t k = 0; k < m_subpop.size(); k++) {
			for (size_t i = 0; i < m_subpop[k].m_first_pop->size(); ++i)
				m_solution[k].push_back(&m_subpop[k].m_first_pop->at(i));
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif // OFEC_DEMO

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::updateDistance(size_t idx, Environment* env) {
		auto pro = env->problem();
		Real min_distance;
		int z = m_source_idx.second;
		for (size_t j = 0; j < m_subpop[idx].m_first_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_subpop.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_CSIWDN(pro)->calculateDistance((*(m_subpop[idx].m_first_pop))[j].variable().index(z), m_subpop[k].m_best->variable().index(z));

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_first_pop))[j].m_distance_fitness = min_distance;
		}
		for (size_t j = 0; j < m_subpop[idx].m_second_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_subpop.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_CSIWDN(pro)->calculateDistance((*(m_subpop[idx].m_second_pop))[j].variable().index(z), m_subpop[k].m_best->variable().index(z));

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_second_pop))[j].m_distance_fitness = min_distance;
		}

		for (size_t j = 0; j < m_subpop[idx].m_first_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_subpop.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_CSIWDN(pro)->calculateDistance((*(m_subpop[idx].m_first_pop))[j].trial().variable().index(z), m_subpop[k].m_best->variable().index(z));

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_first_pop))[j].m_pu_distance_fitness = min_distance;
		}
		for (size_t j = 0; j < m_subpop[idx].m_second_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_subpop.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_CSIWDN(pro)->calculateDistance((*(m_subpop[idx].m_second_pop))[j].trial().variable().index(z), m_subpop[k].m_best->variable().index(z));

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_second_pop))[j].m_pu_distance_fitness = min_distance;
		}
	}

	template <typename TPop1, typename TPop2, typename TInd>
	bool AMP_CC_CSIWDN<TPop1, TPop2, TInd>::isStableByPopNum() {
		if (m_num_be_removed.size() < m_standard_stable) return false;
		int i = m_num_be_removed.size() - 1;
		int count = 0;
		while (m_num_be_removed[i] == 0 && m_num_be_removed[i - 1] == 0) {
			++count;
			--i;
			if (i == 0) break;
		}
		return count >= m_standard_stable - 1;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	bool AMP_CC_CSIWDN<TPop1, TPop2, TInd>::isStableByObj() {
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			if (!m_subpop[i].isStableObj()) {
				return false;
			}
		}
		return true;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	bool AMP_CC_CSIWDN<TPop1, TPop2, TInd>::isStable() {
		if (m_coverage_result.size() < m_standard_stable) return false;
		int q = m_source_idx.first, z = m_source_idx.second;
		if (m_coverage_result.back()[z] == 1) return false;
		int count = 0;
		int i = int(m_coverage_result.size()) - 1;
		bool is_temp_stable = true;
		while (i > 0) {
			for (int k = q; k <= z; ++k) {
				if (m_coverage_result[i][k] != m_coverage_result[i - 1][k]) {
					is_temp_stable = false;
					break;
				}
			}
			if (is_temp_stable)
				++count;
			else
				return count >= m_standard_stable - 1;
			--i;
		}
		return count >= m_standard_stable - 1;
	}


	template <typename TPop1, typename TPop2, typename TInd>
	int AMP_CC_CSIWDN<TPop1, TPop2, TInd>::competition(Environment* env) {
		auto pro = env->problem();
		std::vector<size_t> subpop_indx;
		int q = m_source_idx.first, z = m_source_idx.second;
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			bool is_unique = true;
			if (m_subpop[i].m_best->objective()[0] == 10 )
				continue;
			for (size_t j = 0; j < subpop_indx.size(); ++j) {
				bool is_same_loc = true;
				if (m_subpop[i].m_best->objective()[0]>5) {
					is_unique = false;
					break;
				}
				for (size_t k = q; k <= z; k++) {
					if (m_subpop[subpop_indx[j]].m_best->variable().index(k) != m_subpop[i].m_best->variable().index(k)) {
						is_same_loc = false;
						break;
					}
				}
				if (is_same_loc) {
					m_subpop.populationCompetition(subpop_indx[j], i);
					if (m_subpop[i].m_best->objective()[0] < m_subpop[subpop_indx[j]].m_best->objective()[0])
						*(m_subpop[subpop_indx[j]].m_best) = *(m_subpop[i].m_best);
					is_unique = false;
					break;
				}
			}
			if (is_unique)
				subpop_indx.push_back(i);
		}
		int differ1 = m_subpop.size();
		int differ2 = subpop_indx.size();
		int differ = differ1 - differ2;
		int return_value = differ;
		if (differ > 0) {

			DynSaMultiPop<PopulationType> temp(subpop_indx.size(), m_subpop[0].m_first_pop->size(), env);
			for (size_t i = 0; i < subpop_indx.size(); ++i) {
				temp[i] = m_subpop[subpop_indx[i]];
			}
			for (size_t i = 0; i < subpop_indx.size(); ++i) {
				m_subpop[i] = temp[i];
			}
		}
		while (differ) {
			m_subpop.popBack();
			--differ;
		}
		return return_value;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::addNewPop(size_t num, CSIWDN::InitType type, Environment* env) {
	//	auto pro = env->problem();

		//computeNodeDistance(CAST_CSIWDN(pro)->distanceNode());
		m_subpop.addNewPopulation(env, m_random.get(), num, m_subpop_size, type, m_source_idx);
		m_num_curr_pop += num;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::computeNodeDistance(std::vector<Real> & vect, Environment* env) {
		auto pro = env->problem();
		size_t node_size = CAST_CSIWDN(pro)->numberNode();
		int z = m_source_idx.second, q = m_source_idx.first;
		vect.resize(node_size);
		for (size_t i = 0; i < node_size; ++i) {
			Real curr_best;
			for (size_t j = 0; j < m_subpop.size(); ++j) {
				Real value = CAST_CSIWDN(pro)->calculateDistance(i + 1, m_subpop[j].m_best->variable().index(z));
				if (j == 0) curr_best = value;
				else if (curr_best > value) curr_best = value;
			}
			vect[i] = curr_best;
		}
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void AMP_CC_CSIWDN<TPop1, TPop2, TInd>::output(Environment* env) {
		auto pro = env->problem();

		std::string net = CAST_CSIWDN(pro)->mapName();
		mutex_AMP_out.lock();

		//get time
		time_t t = time(0);
		char tmp[32] = { NULL };
		strftime(tmp, sizeof(tmp), "%d_%H", localtime(&t)); //tmp is string  %Y-%m-%d %H:%M:%S
		//

		std::stringstream path;
		path << g_working_directory << "/result/amp_cc/" << net << "/ahc_best_indi_" << tmp << ".txt";
		std::ofstream fout(path.str(), std::ios::app);
		fout << std::endl;
		fout << "----------------------------------------------------" << std::endl;
		//fout << "runID " << this << ": " << std::endl;
		fout << "maximum evaluations:" << m_maximum_evaluations << std::endl;
		fout << "time window:" << m_eval_time << std::endl;
		fout << "number of subpopulations:" << m_initial_num_pop << std::endl;
		fout << "subpopulation size:" << m_subpop_size << std::endl;
		fout << "number of populations:" << m_subpop.size() << std::endl;
		fout << "use lstm proportion:" << m_lstm_phase << std::endl;
		fout << "----------------------------------------------------" << std::endl;
		fout << std::endl;


		fout << "/***************************************************/" << std::endl;

		std::sort(m_last_sol.begin(), m_last_sol.end(), [](decltype(m_last_sol.front()) a, decltype(m_last_sol.front()) b) {return a.objective(0) < b.objective(0); });
		int size = 3;

		if (m_last_sol.size() < 3)
			size = m_last_sol.size();


		for (size_t i = 0; i < size; i++) {
			for (size_t z = 0; z < CAST_CSIWDN(pro)->numberSource(); z++) {
				fout << "location: " << m_last_sol[i].variable().location(z) << std::endl;
				fout << "start time: " << m_last_sol[i].variable().startTime(z) << std::endl;
				fout << "duration: " << m_last_sol[i].variable().duration(z) << std::endl;
				fout << "multiplier: " << std::endl;
				for (auto& mul_value : m_last_sol[i].variable().multiplier(z))
					fout << mul_value << "  ";
				fout << std::endl;
			}
			fout << std::endl;
			fout << "standard deviation: " << calStandardDeviation(m_last_sol[i].variable(), CAST_CSIWDN(pro)->optima()->solution(0).variable(), env) << std::endl;
			fout << "obj: " << m_last_sol[i].objective()[0] << std::endl;
			fout << std::endl;

			fout << "----------------------------------------------------" << std::endl;

		}
		fout.close();

		std::stringstream path2;
		path2 << g_working_directory << "/result/amp_cc/" << net << "/ahc_curr_best_obj_" << tmp << "h.txt";
		std::ofstream fout2(path2.str(), std::ios::app);
		fout2 << std::endl;
		//fout2 << "runID " << this << ": " << std::endl;
		for (size_t i = 0; i < m_curr_best_obj.size(); ++i)
			fout2 << i << "   " << m_curr_evals[i] << "   " << m_curr_best_obj[i] << std::endl;
		fout2.close();

		std::stringstream path3;
		path3 << g_working_directory << "/result/amp_cc/" << net << "/ahc_coverage_result_" << tmp << "h.txt";
		std::ofstream fout3(path3.str(), std::ios::app);
		fout3 << std::endl;
		//fout3 << "runID " << this << ": " << std::endl;
		for (size_t i = 0; i < m_coverage_result.size(); ++i) {
			fout3 << i << "   " << m_num_pop_result[i] << "   ";
			for (size_t j = 0; j < CAST_CSIWDN(pro)->numberSource(); ++j)
				fout3 << m_coverage_result[i][j] << "  ";
			fout3 << std::endl;
		}
		fout3 << "total_evaluations: " << m_evaluations << std::endl;
		fout3 << "----------------------------------------------------" << std::endl;
		fout3.close();

		std::stringstream path4;
		path4 << g_working_directory << "/result/amp_cc/" << net << "/best_solu_location" << tmp << "h.txt";
		std::ofstream fout4(path4.str(), std::ios::app);
		fout4 << std::endl;
		//fout4 << "runID " << this << ": " << std::endl;
		for (size_t i = 0; i < m_his_best.size(); i++) {
			fout4 << m_curr_evals[i] << " ";
			fout4 << CAST_CSIWDN(pro)->optima()->solution(0).variable().index(m_curr_sourceidx[i]) << " ";
			fout4 << m_his_best[i].variable().index(m_curr_sourceidx[i]) << std::endl;
		}
		fout4 << "----------------------------------------------------" << std::endl;
		fout4.close();

		std::stringstream path5;
		path5 << g_working_directory << "/result/amp_cc/" << net << "/best_solu_injection_rate" << tmp << "h.txt";
		std::ofstream fout5(path5.str(), std::ios::app);
		fout5 << std::endl;
		//fout5 << "runID " << this << ": " << std::endl;
		for (size_t k = 0; k < CAST_CSIWDN(pro)->numberSource(); k++) {
			auto& var = CAST_CSIWDN(pro)->optima()->solution(0).variable();
			for (size_t j = 0; j < var.multiplier(k).size(); ++j)
				fout5 << var.multiplier(k)[j] << " ";
			fout5 << "\n";

			fout5 << m_curr_evals[m_id_his_each_source[k].front()] << " ";
			for (size_t j = 0; j < var.multiplier(k).size(); ++j)
				fout5 << m_his_best[m_id_his_each_source[k].front()].variable().multiplier(k)[j] << " ";
			fout5 << "\n";

			fout5 << m_curr_evals[m_id_his_each_source[k][m_id_his_each_source[k].size() / 6]] << " ";
			for (size_t j = 0; j < var.multiplier(k).size(); ++j)
				fout5 << m_his_best[m_id_his_each_source[k][m_id_his_each_source[k].size() / 6]].variable().multiplier(k)[j] << " ";
			fout5 << "\n";

			fout5 << m_curr_evals[m_id_his_each_source[k][m_id_his_each_source[k].size() / 3]] << " ";
			for (size_t j = 0; j < var.multiplier(k).size(); ++j)
				fout5 << m_his_best[m_id_his_each_source[k][m_id_his_each_source[k].size() / 3]].variable().multiplier(k)[j] << " ";
			fout5 << "\n";

			fout5 << m_curr_evals[m_id_his_each_source[k].back()] << " ";
			for (size_t j = 0; j < var.multiplier(k).size(); ++j)
				fout5 << m_his_best[m_id_his_each_source[k].back()].variable().multiplier(k)[j] << " ";
			fout5 << "\n";

			fout5 << "----------------------------------------------------" << std::endl;
		}
		fout5.close();


		mutex_AMP_out.unlock();
		
	}

}

#endif //! OFEC_AMP_CC_CSIWDN_H