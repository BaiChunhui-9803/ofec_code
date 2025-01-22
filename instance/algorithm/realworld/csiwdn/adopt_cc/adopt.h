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
* L. Liu, S. Ranji Ranjithan, G. Mahinthakumar. Contamination source identification
* in water distribution systems using an adaptive dynamic optimization procedure[J].
* Journal of Water Resources Planning and Management. 2011, 137(2): 183-192.
*********************************************************************************/
// updated Jun 10, 2019 by Li Zhou


#ifndef OFEC_ADOPT_H
#define OFEC_ADOPT_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../dyn_sa_multipop.h"
#include "../cc_framework.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#include<time.h>
#endif // !OFEC_DEMO

namespace ofec {
	template <typename TPop1, typename TPop2, typename TInd>
	class ADOPT : public Algorithm{
		std::mutex mutex_ADOPT_out;
		using PopulationType = CC_framework<TPop1, TPop2, TInd>;

	public:
		void run_();
		void initialize_();
		void initMultiPop();

	private:
		DynSaMultiPop<PopulationType> m_subpop;
		size_t m_subpop_size;
		size_t m_num_curr_pop;
		size_t m_initial_num_pop;
		Real m_tar;
		Real m_tar_multipl;
		int m_iteration = 0;
		size_t m_data_update_eval;
		size_t m_max_iters;
		size_t m_max_eval;
		size_t m_initial_phase;
		size_t m_final_phase;
		long m_eval_time;  //延迟的优化时间
		std::pair<int, int> m_source_idx;
		size_t m_cur_contaminant;
		CSIWDN::ClusterType m_cluster_type;

		std::string m_path_result;
		std::vector<std::vector<Real>> m_coverage_result;
		std::vector<size_t> m_num_pop_result;
		std::vector<Real> m_curr_best_obj;

		std::vector<Solution<VarCSIWDN>> m_last_sol;

		void recordIntermediateResult();
		std::vector<Real> calCoverageRate();
		int evolve();
		size_t elimintePop();
		void addNewPop(size_t num, CSIWDN::InitType type);
		template <typename Population>
		bool isFeasiblePopulation(const Real tar, const Population& pop);
		void updateDistance(size_t idx);
		//bool updateAlternativeSolution(size_t idx, TInd& indi);
		void updateStartEndSouce();
		const std::pair<int, int>& sourceStartEndIdx() { return m_source_idx; }
		void record();
		Real calStandardDeviation(const VarCSIWDN& var, const VarCSIWDN& opt_var);
		void output();


#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	};

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_eval_time = v.get<int>("time window");
		m_num_curr_pop = v.get<int>("number of subpopulations");
		m_initial_num_pop = v.get<int>("number of subpopulations");
		m_subpop_size = v.get<int>("subpopulation size");
		m_tar_multipl = (Real)std::get<double>(v.at("alpha"));
		m_max_iters = v.get<int>("maxIter");
		m_max_eval = v.get<int>("maximum evaluations");
		m_cluster_type = CSIWDN::ClusterType::kNOCluster;
		m_tar = 1.2;
		m_initial_phase = CAST_CSIWDN(m_problem.get())->phase();
		m_final_phase = CAST_CSIWDN(m_problem.get())->totalPhase();
		m_data_update_eval = m_max_eval / (m_final_phase - m_initial_phase + 1);
		m_subpop.clear();
		CAST_CSIWDN(m_problem.get())->setClusterType(m_cluster_type);
		m_cur_contaminant = 0;
		m_source_idx.first = 0;
		m_source_idx.second = 0;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::initMultiPop() {
		CAST_CSIWDN(m_problem.get())->setInitType(CSIWDN::InitType::kRandom);
		CAST_CSIWDN(m_problem.get())->clustering(m_num_curr_pop, m_random.get());
		for (size_t i = 0; i < m_num_curr_pop; ++i) {  // initialize subpopulation 
			auto subpop = std::make_unique<PopulationType>(m_subpop_size, m_problem.get());
			CAST_CSIWDN(m_problem.get())->setPopIdentifier(i);
			subpop->m_first_pop->initialize(m_problem.get(), m_random.get());
			subpop->m_second_pop->initialize(m_problem.get(), m_random.get());
			subpop->fillEach(m_random.get(), m_problem.get(), this, PopulationType::FillType::random, m_source_idx);
			subpop->updateBest(m_problem.get());
			m_subpop.append(subpop);
		}
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::run_() {
		initMultiPop();
#ifdef OFEC_DEMO
		updateBuffer();
#endif // OFEC_DEMO
		int tag = kNormalEval;
		int z = 0;
		CAST_CSIWDN(m_problem.get())->setAlgorithmStart();
		while (!terminating() && tag != kTerminate) {
			for (size_t g = 0; g < m_max_iters; g++)
			{
				tag = evolve();
				recordIntermediateResult();
				m_tar *= (1 - m_tar_multipl);
			}
			m_num_curr_pop = elimintePop();
			//if (m_num_curr_pop == 1 || m_curr_best_obj.back() < 0.00001)
			//	break;
			if (m_evaluations >= (CAST_CSIWDN(m_problem.get())->phase() - m_initial_phase + 1) * m_data_update_eval) {
				int current_phase = m_evaluations / m_data_update_eval + m_initial_phase;
				if (current_phase <= CAST_CSIWDN(m_problem.get())->totalPhase()) {
					CAST_CSIWDN(m_problem.get())->updatePhase(current_phase);
					updateStartEndSouce();
					std::cout << "The current phase:" << current_phase << std::endl;
				}
				for (size_t k = 0; k < m_subpop.size(); ++k) {
					m_subpop[k].evaluate(m_problem.get(), this);
					m_subpop[k].updateBest(m_problem.get());
				}
				if ((z + 1) < CAST_CSIWDN(m_problem.get())->numSource() && CAST_CSIWDN(m_problem.get())->phase() >= ((*(m_subpop[0].m_first_pop))[0].variable().startTime(z + 1) / CAST_CSIWDN(m_problem.get())->intervalTimeStep()) && CAST_CSIWDN(m_problem.get())->phase() < CAST_CSIWDN(m_problem.get())->totalPhase()) {
					addNewPop(m_initial_num_pop, CSIWDN::InitType::kRandom);
					z++;
					std::cout << "The next source" << std::endl;
				}
			}
			if (m_evaluations > m_max_eval)
				break;


			//std::vector<size_t> sub_idx;
			//for (auto& i : m_alternative)
			//	sub_idx.emplace_back(i.second);
			//m_alternative.clear();
			//for (size_t i = 0; i < sub_idx.size(); ++i) {
			//	updateAlternativeSolution(sub_idx[i], *(m_subpop[sub_idx[i]].m_best));
			//}
			//if (m_evaluations >= (CAST_CSIWDN(m_problem.get())->phase() - m_initial_phase + 1) * m_data_update_eval) {
			//	if (m_alternative.size() == 1)
			//		break;
			//	//m_tar *= 0.43;
			//	m_tar *= (1 - m_tar_multipl);
			//	int current_phase = m_evaluations / m_data_update_eval + m_initial_phase;
			//	if (current_phase <= CAST_CSIWDN(m_problem.get())->totalPhase()) {
			//		CAST_CSIWDN(m_problem.get())->updatePhase(current_phase);
			//	}
			//	for (size_t k = 0; k < m_subpop.size(); ++k) {
			//		m_subpop[k].evaluate(m_problem.get(), this);
			//		m_subpop[k].updateBest(m_problem.get());
			//	}
			//	m_alternative.clear();
			//	for (size_t i = 0; i < sub_idx.size(); ++i) {
			//		updateAlternativeSolution(sub_idx[i], *(m_subpop[sub_idx[i]].m_best));
			//	}
			//	if ((z + 1) < CAST_CSIWDN(m_problem.get())->numSource() && CAST_CSIWDN(m_problem.get())->phase() >= ((*(m_subpop[0].m_first_pop))[0].variable().startTime(z + 1) / CAST_CSIWDN(m_problem.get())->intervalTimeStep()) && CAST_CSIWDN(m_problem.get())->phase() < CAST_CSIWDN(m_problem.get())->totalPhase()) {
			//		addNewPop(m_initial_num_pop, CSIWDN::InitType::kRandom);
			//		z++;
			//		std::cout << "The next source" << std::endl;
			//	}
			//}
#ifdef OFEC_DEMO
			updateBuffer();
#endif // OFEC_DEMO
		}
		for (size_t k = 0; k < m_subpop.size(); ++k) {
			m_last_sol.push_back(dynamic_cast<Solution<VarCSIWDN>&>(*m_subpop[k].m_best));
		}
		output();
	}

	template <typename TPop1, typename TPop2, typename TInd>
	size_t ADOPT<TPop1, TPop2, TInd>::elimintePop() {
		std::vector<size_t> subpop_indx;
		int z = m_source_idx.second, q = m_source_idx.first;
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			bool is_unique = true;
			for (size_t j = 0; j < subpop_indx.size(); ++j) {
				bool is_same_loc = true;
				for (size_t k = 0; k <= z; k++) {
					if (m_subpop[subpop_indx[j]].m_best->variable().index(k) != m_subpop[i].m_best->variable().index(k)) {
						is_same_loc = false;
						break;
					}
				}
				if (is_same_loc) {
					if (m_subpop[i].m_best->objective()[0] < m_subpop[subpop_indx[j]].m_best->objective()[0]) {
						subpop_indx[j] = i;
					}
					is_unique = false;
					break;
				}
			}
			if (is_unique)
				subpop_indx.push_back(i);
		}
		size_t differ1 = m_subpop.size();
		size_t differ2 = subpop_indx.size();
		size_t differ = differ1 - differ2;
		if (differ > 0) {
			DynSaMultiPop<PopulationType> temp(subpop_indx.size(), m_subpop[0].m_first_pop->size(), m_problem.get());
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
		return differ2;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::recordIntermediateResult() {
		m_coverage_result.push_back(calCoverageRate());
		m_num_pop_result.push_back(m_subpop.size());
		Real best = m_subpop[0].m_best->objective()[0];
		for (size_t i = 1; i < m_subpop.size(); ++i) {
			if (best > m_subpop[i].m_best->objective()[0])
				best = m_subpop[i].m_best->objective()[0];
		}
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			for (size_t j = 0; j < CAST_CSIWDN(m_problem.get())->numberSource(); ++j) {
				std::cout << m_subpop[i].m_best->variable().location(j) << "  ";
			}
			std::cout << std::endl;
		}
		m_curr_best_obj.push_back(best);
		std::cout << "current best obj:" << best << std::endl;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	std::vector<Real> ADOPT<TPop1, TPop2, TInd>::calCoverageRate() {
		int number_source = CAST_CSIWDN(m_problem.get())->numberSource();
		std::vector <std::vector<bool>> coverage(number_source, std::vector<bool>(CAST_CSIWDN(m_problem.get())->numberNode(), false));
		std::vector<int> count(number_source, 0);
		int z = m_source_idx.second, q = m_source_idx.first;
		for (size_t i = 0; i < m_subpop.size(); ++i) {
			for (size_t j = 0; j < m_subpop[i].m_first_pop->size(); ++j) {
				for (size_t k = q; k <= z; k++) {
					int index = (*(m_subpop[i].m_first_pop))[j].variable().index(k) - 1;
					if (!coverage[k][index]) { coverage[k][index] = true; ++count[k]; }
				}
			}
		}
		std::vector<Real> count_coverage(number_source, 0);
		for (size_t i = q; i <= z; i++) {
			count_coverage[i] = (Real)count[i] / (Real)CAST_CSIWDN(m_problem.get())->numberNode();
		}
		return count_coverage;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	int ADOPT<TPop1, TPop2, TInd>::evolve() {

		int tag = kNormalEval;

		for (size_t i = 0; i < m_subpop.size(); ++i) {
			m_subpop[i].m_first_pop->updateProbability(m_problem.get(), m_source_idx);
			m_subpop[i].m_second_pop->updateCR(m_problem.get(), m_random.get(), m_source_idx);
			m_subpop[i].m_second_pop->updateF(m_problem.get(), m_random.get(), m_source_idx);
			m_subpop[i].m_second_pop->updateProStrategy(m_problem.get(), m_source_idx);

			m_subpop[i].m_first_pop->mutate(m_problem.get(), m_random.get(), m_source_idx);
			m_subpop[i].m_second_pop->mutate(m_problem.get(), m_random.get(), m_source_idx);

			m_subpop[i].m_second_pop->recombine(m_problem.get(), m_random.get(), m_source_idx);


			bool flag_first = isFeasiblePopulation<TPop1>(m_tar, *(m_subpop[i].m_first_pop));
			bool flag_second = isFeasiblePopulation<TPop2>(m_tar, *(m_subpop[i].m_second_pop));
			if (flag_first || flag_second)
				updateDistance(i);
			m_subpop[i].m_first_pop->select(m_problem.get(), this, flag_first, m_source_idx);
			m_subpop[i].m_second_pop->select(m_problem.get(), this, flag_second, m_source_idx);

			m_subpop[i].m_first_pop->iteration()++;
			m_subpop[i].m_second_pop->iteration()++;

			m_subpop[i].fillEach(m_random.get(), m_problem.get(), this, CC_framework<TPop1, TPop2, TInd>::FillType::best, m_source_idx);
			m_subpop[i].updateBest(m_problem.get());
		}
		return tag;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::addNewPop(size_t num, CSIWDN::InitType type) {
		CAST_CSIWDN(m_problem.get())->clustering(m_num_curr_pop + num, m_random.get());
		int q = m_source_idx.first, z = m_source_idx.second;
		for (size_t i = m_num_curr_pop; i < m_num_curr_pop + num; ++i) {
			auto subpop = std::make_unique<PopulationType>(m_subpop_size, m_problem.get());
			CAST_CSIWDN(m_problem.get())->setPopIdentifier(i);
			subpop->m_first_pop->initialize(m_problem.get(), m_random.get());
			subpop->m_second_pop->initialize(m_problem.get(), m_random.get());
			for (size_t j = 0; j < m_subpop_size; j++) {
				size_t I = m_random->uniform.nextNonStd<size_t>(0, m_num_curr_pop);
				for (size_t k = 0; k < z; k++) {
					(*(subpop->m_first_pop))[j].variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
					(*(subpop->m_second_pop))[j].variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
					(*(subpop->m_first_pop))[j].mpu().variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
					(*(subpop->m_second_pop))[j].mpu().variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
					(*(subpop->m_first_pop))[j].mpv().variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
					(*(subpop->m_second_pop))[j].mpv().variable().getEpa(k) = m_subpop[I].m_best->variable().getEpa(k);
				}
			}
			subpop->fillEach(m_random.get(), m_problem.get(), this, PopulationType::FillType::random, m_source_idx);
			subpop->updateBest(m_problem.get());
			m_subpop.append(subpop);
		}
		m_num_curr_pop += num;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	template <typename Population>
	bool ADOPT<TPop1, TPop2, TInd>::isFeasiblePopulation(const Real tar, const Population& pop) {

		size_t phase = CAST_CSIWDN(m_problem.get())->phase();
		size_t interval = CAST_CSIWDN(m_problem.get())->interval();
		size_t num = phase * interval;
		Real temp = 0;

		for (size_t i = 0; i < num; ++i) {
			for (int j = 0; j < CAST_CSIWDN(m_problem.get())->numSensor(); ++j) {
				temp += pow(CAST_CSIWDN(m_problem.get())->observationConcentration()[j][i], 2);
			}
		}

		Real benchmark = tar * sqrt(temp / (CAST_CSIWDN(m_problem.get())->numSensor() * num));

		size_t count_feasible = 0, count_infeasible = 0;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].objective()[0] <= benchmark) ++count_feasible;
			else ++count_infeasible;
		}

		return count_feasible >= count_infeasible;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::updateDistance(size_t idx) {
		Real min_distance, one_dis, dis=0.0;
		int z = m_source_idx.second, q = m_source_idx.first;

		for (size_t i = 0; i < m_subpop[idx].m_first_pop->size(); ++i) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t j = 0; j < m_subpop.size(); ++j) {
				if (idx == j) continue;
				dis = 0.0;
				for (size_t k = q; k <= z; k++) {
					one_dis = CAST_CSIWDN(m_problem.get())->calculateDistance((*(m_subpop[idx].m_first_pop))[i].variable().index(k), m_subpop[j].m_best->variable().index(k));
					dis += one_dis;
				}
				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_first_pop))[i].m_distance_fitness = min_distance;
		}

		for (size_t i = 0; i < m_subpop[idx].m_second_pop->size(); ++i) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t j = 0; j < m_subpop.size(); ++j) {
				if (idx == j) continue;
				dis = 0.0;
				for (size_t k = q; k <= z; k++) {
					one_dis = CAST_CSIWDN(m_problem.get())->calculateDistance((*(m_subpop[idx].m_second_pop))[i].variable().index(k), m_subpop[j].m_best->variable().index(k));
					dis += one_dis;
				}
				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_second_pop))[i].m_distance_fitness = min_distance;
		}

		for (size_t i = 0; i < m_subpop[idx].m_first_pop->size(); ++i) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t j = 0; j < m_subpop.size(); ++j) {
				if (idx == j) continue;
				dis = 0.0;
				for (size_t k = q; k <= z; k++) {
					one_dis = CAST_CSIWDN(m_problem.get())->calculateDistance((*(m_subpop[idx].m_first_pop))[i].trial().variable().index(k), m_subpop[j].m_best->variable().index(k));
					dis += one_dis;
				}
				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_first_pop))[i].m_pu_distance_fitness = min_distance;
		}

		for (size_t i = 0; i < m_subpop[idx].m_second_pop->size(); ++i) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t j = 0; j < m_subpop.size(); ++j) {
				if (idx == j) continue;
				dis = 0.0;
				for (size_t k = q; k <= z; k++) {
					one_dis = CAST_CSIWDN(m_problem.get())->calculateDistance((*(m_subpop[idx].m_second_pop))[i].trial().variable().index(k), m_subpop[j].m_best->variable().index(k));
					dis += one_dis;
				}
				if (min_distance > dis) min_distance = dis;
			}
			(*(m_subpop[idx].m_second_pop))[i].m_pu_distance_fitness = min_distance;
		}

	}

	//template <typename TPop1, typename TPop2, typename TInd>
	//bool ADOPT<TPop1, TPop2, TInd>::updateAlternativeSolution(size_t idx, TInd& indi) {
	//	int z = m_source_idx.second, q = m_source_idx.first;
	//	bool is_same_loc = true;
	//	for (auto& i : m_alternative) {
	//		for (size_t k = 0; k <= z; k++) {
	//			if (i.first->variable().index(k) != indi.variable().index(k)) {
	//				is_same_loc = false;
	//				break;
	//			}
	//		}
	//		if (is_same_loc) {
	//			if (indi.dominate(*(i.first),m_problem.get())) {
	//				i.first.reset(new TInd(indi));
	//				i.second = idx;
	//			}
	//			return false;
	//		}
	//	}
	//	m_alternative.emplace_back(std::make_pair(std::unique_ptr<TInd>(new TInd(indi)), idx));
	//	return true;
	//}

	template<typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::updateStartEndSouce(){
		int q = 0;
		long histry_time = CAST_CSIWDN(m_problem.get())->phase() * CAST_CSIWDN(m_problem.get())->intervalTimeStep() - m_eval_time;
		while (histry_time > 0 && (q + 1) < CAST_CSIWDN(m_problem.get())->numberSource() && histry_time >= (*(m_subpop[0].m_first_pop))[0].variable().startTime(q + 1)) {
			q++;
		}
		m_source_idx.first = q;
		int z = 0;
		while ((z + 1) < CAST_CSIWDN(m_problem.get())->numSource() && CAST_CSIWDN(m_problem.get())->phase() >= ((*(m_subpop[0].m_first_pop))[0].variable().startTime(z + 1) / CAST_CSIWDN(m_problem.get())->intervalTimeStep())) {
			z++;
		}
		m_source_idx.second = z;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::record() {
	}

	template <typename TPop1, typename TPop2, typename TInd>
	Real ADOPT<TPop1, TPop2, TInd>::calStandardDeviation(const VarCSIWDN& var, const VarCSIWDN& opt_var) {
		Real sum = 0;
		size_t total_size = 0;
		for (size_t i = 0; i < CAST_CSIWDN(m_problem.get())->numberSource(); i++) {
			if (var.multiplier(i).size() != opt_var.multiplier(i).size())
				return -1;
			size_t size = var.multiplier(i).size();
			for (size_t j = 0; j < size; j++)
				sum += pow(var.multiplier(i)[j] - opt_var.multiplier(i)[j], 2);
			total_size += size;
		}
		return sqrt(sum / (Real)total_size);
	}

#ifdef OFEC_DEMO
	template <typename TPop1, typename TPop2, typename TInd>
	void ADOPT<TPop1, TPop2, TInd>::updateBuffer() {
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
	void ADOPT<TPop1, TPop2, TInd>::output() {
		std::string net = CAST_CSIWDN(m_problem.get())->mapName();
		mutex_ADOPT_out.lock();
		//get time
		time_t t = time(0);
		char tmp[32] = { NULL };
		strftime(tmp, sizeof(tmp), "%d_%H", localtime(&t)); //tmp is string  %Y-%m-%d %H:%M:%S
		//

		std::stringstream path;
		path << g_working_dir << "/result/adopt/" << net<<"/best_indi_"<<tmp<<".txt";
		auto a = path.str();

		std::ofstream fout(path.str(), std::ios::app);
		fout << std::endl;
		fout << "----------------------------------------------------" << std::endl;
		fout << "runID: " << this << std::endl;
		fout << "maximum evaluations:" << m_max_eval << std::endl;
		fout << "time window:" << m_eval_time << std::endl;
		fout << "number of subpopulations:" << m_initial_num_pop << std::endl;
		fout << "subpopulation size:" << m_subpop_size << std::endl;
		fout << "maxIter:" << m_max_iters << std::endl;
		fout << "----------------------------------------------------" << std::endl;
		fout << std::endl;

		fout << "/***************************************************/" << std::endl;

		std::sort(m_last_sol.begin(), m_last_sol.end(), [](decltype(m_last_sol.front()) a, decltype(m_last_sol.front()) b) {return a.objective(0) < b.objective(0); });
		int size = 3;

		if (m_last_sol.size() < 3)
			size = m_last_sol.size();

		for (size_t i = 0; i < size; ++i) {
			fout << "best solution " << i + 1 << ": " << std::endl;
			for (size_t k = 0; k < CAST_CSIWDN(m_problem.get())->numberSource(); k++) {
				fout << "location: " << m_last_sol[i].variable().location(k) << std::endl;
				fout << "start time: " << m_last_sol[i].variable().startTime(k) << std::endl;
				fout << "duration: " << m_last_sol[i].variable().duration(k) << std::endl;
				fout << "multiplier: " << std::endl;
				for (auto& mul_value : m_last_sol[i].variable().multiplier(k))
					fout << mul_value << "  ";
				fout << std::endl;
			}
			fout << std::endl;
			fout << "standard deviation: " << calStandardDeviation(m_last_sol[i].variable(), CAST_CSIWDN(m_problem.get())->optima()->variable(0)) << std::endl;
			fout << "obj: " << m_last_sol[i].objective()[0] << std::endl;
			fout << std::endl;

			fout << "----------------------------------------------------" << std::endl;

		}
		fout.close();

		std::stringstream path2;
		path2 << g_working_dir << "/result/adopt/" << net << "/coverage_result_" << tmp << ".txt";
		std::ofstream fout2(path2.str(), std::ios::app);
		fout2 << std::endl;
		fout2 << "runID " << this << ": " << std::endl;
		for (size_t i = 0; i < m_coverage_result.size(); ++i) {
			fout2 << i << "   " << m_num_pop_result[i] << "   ";
			for (size_t j = 0; j < CAST_CSIWDN(m_problem.get())->numberSource(); ++j)
				fout2 << m_coverage_result[i][j] << "  ";
			fout2 << std::endl;
		}
		fout2 << "total_evaluations: " << m_evaluations << std::endl;
		fout2 << "----------------------------------------------------" << std::endl;
		fout2.close();



		std::stringstream path3;
		path3 << g_working_dir << "/result/adopt/" << net << "/curr_best_obj_" << tmp << ".txt";
		std::ofstream fout3(path3.str(), std::ios::app);
		fout3 << "runID " << this << ": " << std::endl;
		for (size_t i = 0; i < m_curr_best_obj.size(); ++i)
			fout3 << i << "   " << m_curr_best_obj[i] << std::endl;
		fout3.close();

		mutex_ADOPT_out.unlock();
	}



	//原代码

	////std::mutex mutex_ADOPT_out;
	////template <typename TPop1, typename TPop2, typename TInd>
	////ADOPT<TPop1, TPop2, TInd>::ADOPT(param_map & v) : algorithm(v.at("algorithm name")), \
	////	m_subpop((size_t)v.at("numSubPop"), (size_t)v.at("subpopulation size")), m_tar((Real)1.2), m_tar_multiplier((Real)v.at("alpha")), \
	////	m_path_result(v.at("dataFileResult"))
	////{
	////	m_max_eval = v.at("maximum evaluations");
	////	m_initial_phase = CAST_RP_EPANET->phase();
	////	m_final_phase = CAST_RP_EPANET->m_total_phase;
	////	m_data_update_eval = m_max_eval / (m_final_phase - m_initial_phase + 1);
	////}

	////template <typename TPop1, typename TPop2, typename TInd>
	////void ADOPT<TPop1, TPop2, TInd>::initialize() {
	////	for (size_t i = 0; i < m_subpop.size(); ++i) {  // initialize subpopulation 
	////		CAST_RP_EPANET->set_init_type(init_type::random);
	////		m_subpop[i].m_first_pop->initialize();
	////		m_subpop[i].m_second_pop->initialize();
	////		m_subpop[i].fill_each(CC_framework<TPop1, TPop2, TInd>::fill_type::random);
	////		m_subpop[i].update_best();
	////		m_alternative.emplace_back(std::make_pair(std::unique_ptr<TInd>(new TInd(*(m_subpop[i].m_best))), i));
	////	}
	////}

	//template <typename TPop1, typename TPop2, typename TInd>
	//void ADOPT<TPop1, TPop2, TInd>::run_() {
	//	int tag = kNormalEval;

	//	while (tag != kTerminate) {

	//		//size_t total_evaluations = CAST_RP_EPANET->total_evaluations();

	//		record_intermediate_result();

	//		//std::cout << "runID: " << global::ms_global.get()->m_runID << " iteraton: " << m_iteration << "  phase: " << CAST_RP_EPANET->phase() << "  total_eval: " << total_evaluations << "  alternative: " << m_alternative.size() << " coverage: " << m_coverage_result[m_iteration] << std::endl;

	//		tag = evolve();

	//		size_t curr_eval = global::ms_global->m_problem->evaluations();
	//		if (curr_eval >= (CAST_RP_EPANET->phase() - m_initial_phase + 1) * m_data_update_eval) {

	//			std::vector<size_t> sub_idx;
	//			for (auto &i : m_alternative)
	//				sub_idx.emplace_back(i.second);
	//			m_alternative.clear();
	//			for (size_t i = 0; i < sub_idx.size(); ++i) {
	//				update_alternative_solution(sub_idx[i], *(m_subpop[sub_idx[i]].m_best));
	//			}

	//			if (m_alternative.size() == 1) break;
	//			//m_tar *= 0.43;
	//			m_tar *= m_tar_multiplier;
	//			CAST_RP_EPANET->phase_next();
	//			if (CAST_RP_EPANET->is_time_out()) break;
	//		}

	//	}
	//	output();
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//EvaluationTag ADOPT<TPop1, TPop2, TInd>::evolve() {

	//	int tag = kNormalEval;
	//	int count = 0;

	//	for (size_t i = 0; i < m_alternative.size(); ++i) {
	//		size_t idx = m_alternative[i].second;

	//		m_subpop[idx].m_first_pop->update_probability();
	//		m_subpop[idx].m_second_pop->update_CR();
	//		m_subpop[idx].m_second_pop->update_F();
	//		m_subpop[idx].m_second_pop->update_pro_strategy();

	//		m_subpop[idx].m_first_pop->mutate();
	//		m_subpop[idx].m_second_pop->mutate();

	//		m_subpop[idx].m_second_pop->recombine();

	//		bool flag_first = is_feasible_population<TPop1>(m_tar, *(m_subpop[idx].m_first_pop));
	//		bool flag_second = is_feasible_population<TPop2>(m_tar, *(m_subpop[idx].m_second_pop));
	//		if (flag_first || flag_second)
	//			update_distance(idx);
	//		m_subpop[idx].m_first_pop->select(flag_first);
	//		m_subpop[idx].m_second_pop->select(flag_second);

	//		m_subpop[idx].fillEach(CC_framework<TPop1, TPop2, TInd>::fill_type::best);
	//	}
	//	++m_iteration;
	//	return tag;
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//void ADOPT<TPop1, TPop2, TInd>::update_distance(size_t idx) {
	//	Real min_distance;

	//	for (size_t j = 0; j < m_subpop[idx].m_first_pop->size(); ++j) {
	//		min_distance = std::numeric_limits<Real>::max();
	//		for (size_t k = 0; k < m_subpop.size(); ++k) {
	//			if (idx == k) continue;
	//			Real dis = CAST_RP_EPANET->calculate_distance((*(m_subpop[idx].m_first_pop))[j].variable().index(), m_subpop[k].m_best->variable().index());

	//			if (min_distance > dis) min_distance = dis;
	//		}
	//		(*(m_subpop[idx].m_first_pop))[j].m_distance_fitness = min_distance;
	//	}
	//	for (size_t j = 0; j < m_subpop[idx].m_second_pop->size(); ++j) {
	//		min_distance = std::numeric_limits<Real>::max();
	//		for (size_t k = 0; k < m_subpop.size(); ++k) {
	//			if (idx == k) continue;
	//			Real dis = CAST_RP_EPANET->calculate_distance((*(m_subpop[idx].m_second_pop))[j].variable().index(), m_subpop[k].m_best->variable().index());

	//			if (min_distance > dis) min_distance = dis;
	//		}
	//		(*(m_subpop[idx].m_second_pop))[j].m_distance_fitness = min_distance;
	//	}

	//	for (size_t j = 0; j < m_subpop[idx].m_first_pop->size(); ++j) {
	//		min_distance = std::numeric_limits<Real>::max();
	//		for (size_t k = 0; k < m_subpop.size(); ++k) {
	//			if (idx == k) continue;
	//			Real dis = CAST_RP_EPANET->calculate_distance((*(m_subpop[idx].m_first_pop))[j].trial().variable().index(), m_subpop[k].m_best->variable().index());

	//			if (min_distance > dis) min_distance = dis;
	//		}
	//		(*(m_subpop[idx].m_first_pop))[j].m_pu_distance_fitness = min_distance;
	//	}
	//	for (size_t j = 0; j < m_subpop[idx].m_second_pop->size(); ++j) {
	//		min_distance = std::numeric_limits<Real>::max();
	//		for (size_t k = 0; k < m_subpop.size(); ++k) {
	//			if (idx == k) continue;
	//			Real dis = CAST_RP_EPANET->calculate_distance((*(m_subpop[idx].m_second_pop))[j].trial().variable().index(), m_subpop[k].m_best->variable().index());

	//			if (min_distance > dis) min_distance = dis;
	//		}
	//		(*(m_subpop[idx].m_second_pop))[j].m_pu_distance_fitness = min_distance;
	//	}

	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//bool ADOPT<TPop1, TPop2, TInd>::update_alternative_solution(size_t idx, TInd & indi) {
	//	for (auto& i : m_alternative) {
	//		if (i.first->variable().index() == indi.variable().index()) {
	//			if (indi.dominate(*(i.first))) {
	//				//i.first.release();
	//				i.first.reset(new TInd(indi));
	//				i.second = idx;
	//			}
	//			return false;
	//		}
	//	}
	//	m_alternative.emplace_back(std::make_pair(std::unique_ptr<TInd>(new TInd(indi)), idx));
	//	return true;
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//void ADOPT<TPop1, TPop2, TInd>::record() {
	//	
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//Real ADOPT<TPop1, TPop2, TInd>::cal_standard_deviation(const std::vector<float>& vec1, const std::vector<float>& vec2) {
	//	if (vec1.size() != vec2.size()) return -1;
	//	int size = vec1.size();
	//	Real sum = 0;
	//	for (size_t i = 0; i < size; ++i) {
	//		sum += pow(vec1[i] - vec2[i], 2);
	//	}
	//	return sqrt(sum / (Real)size);
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//Real ADOPT<TPop1, TPop2, TInd>::cal_coverage_rate() {
	//	std::vector<bool> coverage(CAST_RP_EPANET->number_node(), false);
	//	int count = 0;
	//	for (size_t i = 0; i < m_alternative.size(); ++i) {
	//		int idx = m_alternative[i].second;
	//		for (size_t j = 0; j < m_subpop[idx].m_first_pop->size(); ++j) {
	//			int index = (*(m_subpop[idx].m_first_pop))[j].variable().index() - 1;
	//			if (!coverage[index]) { coverage[index] = true; ++count; }
	//		}
	//	}
	//	return (Real)count / (Real)CAST_RP_EPANET->number_node();
	//}

	////template <typename TPop1, typename TPop2, typename TInd>
	////void ADOPT<TPop1, TPop2, TInd>::record_intermediate_result() {
	////	m_coverage_result.push_back(cal_coverage_rate());
	////	m_num_pop_result.push_back(m_alternative.size());
	////	Real best = m_alternative[0].first->objective()[0];
	////	for (size_t i = 1; i < m_alternative.size(); ++i) {
	////		if (best > m_alternative[i].first->objective()[0])
	////			best = m_alternative[i].first->objective()[0];
	////	}
	////	m_curr_best_obj.push_back(best);
	////}

	//template <typename TPop1, typename TPop2, typename TInd>
	//void ADOPT<TPop1, TPop2, TInd>::output() {
	//	mutex_ADOPT_out.lock();
	//	std::stringstream path;
	//	//path << "F:/A_Li_Zhou/Experiment/污水源追踪问题/NET6/40/ADOPT_CCSaDE_NoDistEva_1203.txt";
	//	path << g_working_dir << "/result/new/Alg4Epa/ADOPT/" << m_path_result;

	//	std::ofstream fout(path.str(), std::ios::app);
	//	fout << std::endl;

	//	fout << "runID " << global::ms_global.get()->m_runID << ": " << std::endl;

	//	for (size_t i = 0; i < m_alternative.size(); ++i) {
	//		fout << "best solution " << i + 1 << ": " << std::endl;
	//		fout << "location: " << m_alternative[i].first->variable().location() << std::endl;
	//		fout << "start time: " << m_alternative[i].first->variable().start_time() << std::endl;
	//		fout << "duration: " << m_alternative[i].first->variable().duration() << std::endl;
	//		fout << "multiplier: " << std::endl;
	//		for (auto &mul_value : m_alternative[i].first->variable().multiplier())
	//			fout << mul_value << "  ";
	//		fout << std::endl;
	//		fout << "obj: " << m_alternative[i].first->objective()[0] << std::endl;
	//		fout << "standard deviation: " << cal_standard_deviation(m_alternative[i].first->variable().multiplier(), CAST_RP_EPANET->get_optima().variable(0).multiplier()) << std::endl;
	//		fout << std::endl;
	//		fout << "concentration: " << std::endl;
	//		std::vector<std::vector<float>> concen_data(4, std::vector<float>(144));
	//		CAST_RP_EPANET->get_data(m_alternative[i].first->variable(), concen_data);
	//		for (size_t j = 0; j < 144; ++j) {
	//			fout << concen_data[0][j] << " " << concen_data[1][j] << " " << concen_data[2][j] << " " << concen_data[3][j] << " " << (concen_data[0][j] + concen_data[1][j] + concen_data[2][j] + concen_data[3][j]) / 4.0 << std::endl;
	//		}
	//	}

	//	for (size_t i = 0; i<m_coverage_result.size(); ++i)
	//		fout << i << "   " << m_num_pop_result[i] << "   " << m_coverage_result[i] << std::endl;
	//	int total_evaluations = CAST_RP_EPANET->total_evaluations();
	//	fout << "total_evaluations: " << total_evaluations << std::endl;
	//	fout << "----------------------------------------------------" << std::endl;
	//	fout.close();

	//	std::stringstream path2;
	//	path2 << g_working_dir << "/result/new/Alg4Epa/ADOPT/NET2_97/20w/curr_best_obj.txt";
	//	std::ofstream fout2(path2.str(), std::ios::app);
	//	for (size_t i = 0; i<m_curr_best_obj.size(); ++i)
	//		fout2 << i << "   " << m_curr_best_obj[i] << std::endl;
	//	fout2.close();

	//	mutex_ADOPT_out.unlock();
	//}

	//template <typename TPop1, typename TPop2, typename TInd>
	//template <typename Population>
	//bool ADOPT<TPop1, TPop2, TInd>::is_feasible_population(const Real tar, const Population & pop) {

	//	size_t phase = CAST_RP_EPANET->phase();
	//	size_t interval = CAST_RP_EPANET->interval();
	//	size_t num = phase*interval;
	//	Real temp = 0;

	//	for (size_t i = 0; i < num; ++i) {
	//		for (int j = 0; j < CAST_RP_EPANET->num_sensor(); ++j) {
	//			temp += pow(CAST_RP_EPANET->observation_concentration()[j][i], 2);
	//		}
	//	}

	//	Real benchmark = tar * sqrt(temp / (CAST_RP_EPANET->num_sensor()*num));

	//	size_t count_feasible = 0, count_infeasible = 0;
	//	for (size_t i = 0; i < pop.size(); ++i) {
	//		if (pop[i].objective()[0] <= benchmark) ++count_feasible;
	//		else ++count_infeasible;
	//	}

	//	return count_feasible >= count_infeasible;
	//}

}
#endif // OFEC_ADOPT_CCSaDE_H

