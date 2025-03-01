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
* Li Liu,Emily M. Zechman, G. Mahinthakumar & S. Ranji Ranjithan. Identifying
* contaminant sources for water distribution systems using a hybrid method[J].
* Journal of Water Resources Planning and Management. 2011, 137(2): 183-192.
*********************************************************************************/
// updated Jun 10, 2019 by Li Zhou


#ifndef OFEC_LRM_ADOPT_H
#define OFEC_LRM_ADOPT_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../DynSaMultiPop.h"
#include "../CC_framework.h"
#include "../../../../../../utility/logistic_regression.h"

namespace ofec {
	template <typename Population1, typename Population2, typename Solution>
	class LRM_ADOPT : public algorithm
	{
		std::mutex mutex_LRM_ADOPT_out;
		//using fill_type = CC_framework<Population1, Population2, Solution>::fill_type;
		using PAIR = std::pair<int, Real>;
		struct CmpByValue {
			bool operator()(const PAIR& lhs, const PAIR& rhs) {
				return lhs.second < rhs.second;
			}
		};
	public:
		LRM_ADOPT(param_map & v);
		void run_();
		void initialize();
#ifdef OFEC_DEMO
		void updateBuffer() override{}
#endif
	protected:
		int evolve();

		void update_distance(size_t idx);
		bool update_alternative_solution(size_t idx, Solution & indi);
		template <typename Population>
		bool is_feasible_population(const Real tar, const Population & cc);
		Real cal_standard_deviation(const std::vector<float>& vec1, const std::vector<float>& vec2);
		Real cal_coverage_rate();
		void record_intermediate_result();
		void generate_training_data();
		void establish_logistic_regression(Real  rate, int iter);   // LRM
		void update_index_flag();
		template <typename Population>
		void local_search(Population & pop);  // LS

		std::vector<Real> D2_transform_D1(const std::vector<std::vector<Real>> data);
		void output();
		void record();
	private:
		Real m_tar;
		DynSaMultiPop<CC_framework<Population1, Population2, Solution>> m_sub_population;
		std::vector<std::pair<std::unique_ptr<Solution>, size_t>> m_alternative;
		size_t m_phase_interval;
		Real m_tar_multiplier;
		int m_iteration = 0;
		const std::string m_path_result;

		std::vector<Real> m_coverage_result;
		std::vector<size_t> m_num_pop_result;
		std::vector<Real> m_curr_best_obj;

		std::vector<std::vector<std::vector<Real>>> m_training_data;
		std::vector<std::vector<Real>> m_pro_index_is_real;

		Real m_learning_rate;
		int m_iter_LR;

		size_t m_data_update_eval;
		size_t m_max_eval;
		size_t m_initial_phase;
		size_t m_final_phase;
	};

	//std::mutex mutex_LRM_ADOPT_out;

	template <typename Population1, typename Population2, typename Solution>
	LRM_ADOPT<Population1, Population2, Solution>::LRM_ADOPT(param_map & v) : algorithm(v.at("algorithm name")), \
		m_sub_population((size_t)v.at("numSubPop"), (size_t)v.at("subpopulation size")), m_tar((Real)1.2), m_tar_multiplier((Real)v.at("alpha")), \
		m_path_result(v.at("dataFileResult")), m_learning_rate(v.at("case")), m_iter_LR(v.at("gamma"))
	{
		m_max_eval = v.at("maximum evaluations");
		m_initial_phase = CAST_RP_EPANET->phase();
		m_final_phase = CAST_RP_EPANET->m_total_phase;
		m_data_update_eval = m_max_eval / (m_final_phase - m_initial_phase + 1);
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::initialize() {
		generate_training_data();
		establish_logistic_regression(m_learning_rate, m_iter_LR);
		update_index_flag();
		for (size_t i = 0; i < m_sub_population.size(); ++i) {  // initialize subpopulation 

			CAST_RP_EPANET->set_init_type(init_type::random);

			m_sub_population[i].m_first_pop->initialize();

			m_sub_population[i].m_second_pop->initialize();

			m_sub_population[i].fill_each(CC_framework<Population1, Population2, Solution>::fill_type::random);

			m_sub_population[i].update_best();

			m_alternative.emplace_back(std::make_pair(std::unique_ptr<Solution>(new Solution(*(m_sub_population[i].m_best))), i));
		}
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::run_() {
		int tag = kNormalEval;

		while (tag != kTerminate) {

			//size_t total_evaluations = CAST_RP_EPANET->total_evaluations();

			record_intermediate_result();

			//std::cout << "runID: " << global::ms_global.get()->m_runID << " iteraton: " << m_iteration << "  phase: " << CAST_RP_EPANET->phase() << "  total_eval: " << total_evaluations << "  alternative: " << m_alternative.size() << " coverage: " << m_coverage_result[m_iteration] << std::endl;

			update_index_flag();

			tag = evolve();

			size_t curr_eval = global::ms_global->m_problem->evaluations();
			if (curr_eval >= (CAST_RP_EPANET->phase() - m_initial_phase + 1) * m_data_update_eval) {

				std::vector<size_t> sub_idx;
				for (auto &i : m_alternative)
					sub_idx.emplace_back(i.second);
				m_alternative.clear();
				for (size_t i = 0; i < sub_idx.size(); ++i) {
					update_alternative_solution(sub_idx[i], *(m_sub_population[sub_idx[i]].m_best));
				}

				if (m_alternative.size() == 1) break;
				//m_tar *= 0.43;
				m_tar *= m_tar_multiplier;
				CAST_RP_EPANET->phase_next();
				if (CAST_RP_EPANET->is_time_out()) break;
			}

		}
		output();
	}

	template <typename Population1, typename Population2, typename Solution>
	EvaluationTag LRM_ADOPT<Population1, Population2, Solution>::evolve() {

		int tag = kNormalEval;
		int count = 0;

		for (size_t i = 0; i < m_alternative.size(); ++i) {
			size_t idx = m_alternative[i].second;

			m_sub_population[idx].m_first_pop->update_probability();
			//m_sub_population[idx].m_second_pop->update_probability();
			m_sub_population[idx].m_second_pop->update_CR();
			m_sub_population[idx].m_second_pop->update_F();
			m_sub_population[idx].m_second_pop->update_pro_strategy();

			m_sub_population[idx].m_first_pop->mutate();
			m_sub_population[idx].m_second_pop->mutate();

			m_sub_population[idx].m_second_pop->recombine();

			/*if (is_feasible_population<GL_population>(m_tar, *(m_sub_population[idx].m_first_pop))) {
			update_distance(idx);
			m_sub_population[idx].m_first_pop->select(true);
			m_sub_population[idx].m_second_pop->select(true);
			}
			else {
			m_sub_population[idx].m_first_pop->select(false);
			m_sub_population[idx].m_second_pop->select(false);
			}*/

			/*m_sub_population[i].m_first_pop->select(false);
			m_sub_population[i].m_second_pop->select(false);*/

			bool flag_first = is_feasible_population<Population1>(m_tar, *(m_sub_population[idx].m_first_pop));
			bool flag_second = is_feasible_population<Population2>(m_tar, *(m_sub_population[idx].m_second_pop));
			if (flag_first || flag_second)
				update_distance(idx);
			m_sub_population[idx].m_first_pop->select(flag_first);
			m_sub_population[idx].m_second_pop->select(flag_second);

			local_search(*(m_sub_population[idx].m_first_pop));
			local_search(*(m_sub_population[idx].m_second_pop));

			m_sub_population[idx].fill_each(CC_framework<Population1, Population2, Solution>::fill_type::best);
		}
		++m_iteration;
		return tag;
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::update_distance(size_t idx) {
		Real min_distance;

		for (size_t j = 0; j < m_sub_population[idx].m_first_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_sub_population.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_RP_EPANET->calculate_distance((*(m_sub_population[idx].m_first_pop))[j].variable().index(), m_sub_population[k].m_best->variable().index());

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_sub_population[idx].m_first_pop))[j].m_distance_fitness = min_distance;
		}
		for (size_t j = 0; j < m_sub_population[idx].m_second_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_sub_population.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_RP_EPANET->calculate_distance((*(m_sub_population[idx].m_second_pop))[j].variable().index(), m_sub_population[k].m_best->variable().index());

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_sub_population[idx].m_second_pop))[j].m_distance_fitness = min_distance;
		}

		for (size_t j = 0; j < m_sub_population[idx].m_first_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_sub_population.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_RP_EPANET->calculate_distance((*(m_sub_population[idx].m_first_pop))[j].trial().variable().index(), m_sub_population[k].m_best->variable().index());

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_sub_population[idx].m_first_pop))[j].m_pu_distance_fitness = min_distance;
		}
		for (size_t j = 0; j < m_sub_population[idx].m_second_pop->size(); ++j) {
			min_distance = std::numeric_limits<Real>::max();
			for (size_t k = 0; k < m_sub_population.size(); ++k) {
				if (idx == k) continue;
				Real dis = CAST_RP_EPANET->calculate_distance((*(m_sub_population[idx].m_second_pop))[j].trial().variable().index(), m_sub_population[k].m_best->variable().index());

				if (min_distance > dis) min_distance = dis;
			}
			(*(m_sub_population[idx].m_second_pop))[j].m_pu_distance_fitness = min_distance;
		}

	}

	template <typename Population1, typename Population2, typename Solution>
	bool LRM_ADOPT<Population1, Population2, Solution>::update_alternative_solution(size_t idx, Solution & indi) {
		for (auto& i : m_alternative) {
			if (i.first->variable().index() == indi.variable().index()) {
				if (indi.dominate(*(i.first))) {
					//i.first.release();
					i.first.reset(new Solution(indi));
					i.second = idx;
				}
				return false;
			}
		}
		m_alternative.emplace_back(std::make_pair(std::unique_ptr<Solution>(new Solution(indi)), idx));
		return true;
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::record() {
		//measure::get_measure()->record(global::ms_global.get(), 0, 0);
	}

	template <typename Population1, typename Population2, typename Solution>
	Real LRM_ADOPT<Population1, Population2, Solution>::cal_standard_deviation(const std::vector<float>& vec1, const std::vector<float>& vec2) {
		if (vec1.size() != vec2.size()) return -1;
		int size = vec1.size();
		Real sum = 0;
		for (size_t i = 0; i < size; ++i) {
			sum += pow(vec1[i] - vec2[i], 2);
		}
		return sqrt(sum / (Real)size);
	}

	template <typename Population1, typename Population2, typename Solution>
	Real LRM_ADOPT<Population1, Population2, Solution>::cal_coverage_rate() {
		std::vector<bool> coverage(CAST_RP_EPANET->number_node(), false);
		int count = 0;
		for (size_t i = 0; i < m_alternative.size(); ++i) {
			int idx = m_alternative[i].second;
			for (size_t j = 0; j < m_sub_population[idx].m_first_pop->size(); ++j) {
				int index = (*(m_sub_population[idx].m_first_pop))[j].variable().index() - 1;
				if (!coverage[index]) { coverage[index] = true; ++count; }
			}
		}
		return (Real)count / (Real)CAST_RP_EPANET->number_node();
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::record_intermediate_result() {
		m_coverage_result.push_back(cal_coverage_rate());
		m_num_pop_result.push_back(m_alternative.size());
		Real best = m_alternative[0].first->objective()[0];
		for (size_t i = 1; i < m_alternative.size(); ++i) {
			if (best > m_alternative[i].first->objective()[0])
				best = m_alternative[i].first->objective()[0];
		}
		m_curr_best_obj.push_back(best);
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::output() {
		mutex_LRM_ADOPT_out.lock();
		std::stringstream path;
		path << g_working_dir << "/result/new/Alg4Epa/LRM_ADOPT/" << m_path_result;

		std::ofstream fout(path.str(), std::ios::app);
		fout << std::endl;

		fout << "runID " << global::ms_global.get()->m_runID << ": " << std::endl;

		for (size_t i = 0; i < m_alternative.size(); ++i) {
			fout << "best solution " << i + 1 << ": " << std::endl;
			fout << "location: " << m_alternative[i].first->variable().location() << std::endl;
			fout << "start time: " << m_alternative[i].first->variable().start_time() << std::endl;
			fout << "duration: " << m_alternative[i].first->variable().duration() << std::endl;
			fout << "multiplier: " << std::endl;
			for (auto &mul_value : m_alternative[i].first->variable().multiplier())
				fout << mul_value << "  ";
			fout << std::endl;
			fout << "obj: " << m_alternative[i].first->objective()[0] << std::endl;
			fout << "standard deviation: " << cal_standard_deviation(m_alternative[i].first->variable().multiplier(), CAST_RP_EPANET->get_optima().variable(0).multiplier()) << std::endl;
			fout << std::endl;
			fout << "concentration: " << std::endl;
			std::vector<std::vector<float>> concen_data(4, std::vector<float>(144));
			CAST_RP_EPANET->get_data(m_alternative[i].first->variable(), concen_data);
			for (size_t j = 0; j < 144; ++j) {
				fout << concen_data[0][j] << " " << concen_data[1][j] << " " << concen_data[2][j] << " " << concen_data[3][j] << " " << (concen_data[0][j] + concen_data[1][j] + concen_data[2][j] + concen_data[3][j]) / 4.0 << std::endl;
			}
		}

		for (size_t i = 0; i<m_coverage_result.size(); ++i)
			fout << i << "   " << m_num_pop_result[i] << "   " << m_coverage_result[i] << std::endl;
		int total_evaluations = CAST_RP_EPANET->total_evaluations();
		fout << "total_evaluations: " << total_evaluations << std::endl;
		fout << "----------------------------------------------------" << std::endl;
		fout.close();

		std::stringstream path2;
		path2 << g_working_dir << "/result/new/Alg4Epa/LRM_ADOPT/NET2_97/20w/curr_best_obj.txt";
		std::ofstream fout2(path2.str(), std::ios::app);
		for (size_t i = 0; i<m_curr_best_obj.size(); ++i)
			fout2 << i << "   " << m_curr_best_obj[i] << std::endl;
		fout2.close();

		mutex_LRM_ADOPT_out.unlock();
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::generate_training_data() {

		int start = CAST_RP_EPANET->m_first_phase;
		int end = CAST_RP_EPANET->m_total_phase;
		int num_train = CAST_RP_EPANET->number_node();
		int num_sensor = CAST_RP_EPANET->num_sensor();
		long step = CAST_RP_EPANET->pattern_step();
		int interval = CAST_RP_EPANET->interval();
		int total_unit_time = CAST_RP_EPANET->total_run_time() / CAST_RP_EPANET->m_quality_time_step;

		CAST_RP_EPANET->set_phase(end);

		m_training_data.resize(num_train, std::vector<std::vector<Real>>(end - start + 1, std::vector<Real>(num_sensor))); // node<time<sensor<Real>>>
		variable_epanet sol;
		for (size_t i = 0; i < num_train; ++i) {
			sol.index() = i + 1;
			sol.flag_location() = false;
			sol.start_time() = step * global::ms_global->m_uniform[caller::Problem]->next_non_standard<long>
				(CAST_RP_EPANET->m_min_start_time / step, CAST_RP_EPANET->m_max_start_time / step);
			sol.duration() = const_cast<variable_epanet &>(CAST_RP_EPANET->get_optima().variable(0)).duration();
			size_t size = sol.duration() / step;
			if (sol.duration() % step != 0) ++size;
			sol.multiplier().resize(size);
			for (auto &i : sol.multiplier()) {
				i = global::ms_global->m_uniform[caller::Problem]->next_non_standard<float>(CAST_RP_EPANET->m_min_multiplier, CAST_RP_EPANET->m_max_multiplier);
			}

			std::vector<std::vector<float>> temp_data(num_sensor, std::vector<float>(total_unit_time));
			CAST_RP_EPANET->get_data(sol, temp_data);

			int count = -1;
			for (size_t j = (start - 1) * interval; j <= (end - 1) * interval; j += interval) {
				++count;
				for (size_t k = 0; k < num_sensor; ++k) {
					m_training_data[i][count][k] = temp_data[k][j];
				}
			}
		}

		CAST_RP_EPANET->set_phase(start);
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::establish_logistic_regression(Real  rate, int iter) {

		int num_node = CAST_RP_EPANET->number_node();
		int num_sensor = CAST_RP_EPANET->num_sensor();
		int index_real = CAST_RP_EPANET->get_optima().variable(0).index();
		std::vector<Real> labels(num_node, 0);
		labels[index_real - 1] = 1;

		int size_time = CAST_RP_EPANET->m_total_phase - CAST_RP_EPANET->m_first_phase + 1;

		m_pro_index_is_real.resize(size_time, std::vector<Real>(num_node, 0));

		for (size_t i = 0; i < size_time; ++i) {
			logistic_regression<Real> lr;
			std::vector<std::vector<Real>> training_data_2D(m_training_data[i]);
			std::vector<Real> training_data_1D(D2_transform_D1(training_data_2D));
			lr.init(training_data_1D.data(), labels.data(), num_node, num_sensor, -1, rate, iter, 0, 1);
			lr.train();
			for (size_t j = 0; j < num_node; ++j) {
				m_pro_index_is_real[i][j] = lr.predict(m_training_data[j][i].data(), num_sensor);
			}
		}

	}

	template <typename Population1, typename Population2, typename Solution>
	std::vector<Real> LRM_ADOPT<Population1, Population2, Solution>::D2_transform_D1(const std::vector<std::vector<Real>> data) {
		int size1 = data.size();
		int size2 = data[0].size();
		std::vector<Real> temp(size1*size2);
		for (size_t i = 0; i < size1; ++i) {
			for (size_t j = 0; j < size2; ++j) {
				temp.push_back(data[i][j]);
			}
		}
		return temp;
	}

	template <typename Population1, typename Population2, typename Solution>
	void LRM_ADOPT<Population1, Population2, Solution>::update_index_flag() {
		int num_node = CAST_RP_EPANET->number_node();
		int time_start = CAST_RP_EPANET->m_first_phase - 1;
		int now_phase = CAST_RP_EPANET->phase();

		CAST_RP_EPANET->m_index_flag.resize(num_node, true);

		for (size_t i = 0; i < num_node; ++i) {
			if (m_pro_index_is_real[now_phase - time_start - 1][i] == 0)
				CAST_RP_EPANET->m_index_flag[i] = false;
		}
	}

	template <typename Population1, typename Population2, typename Solution>
	template <typename Population>
	bool LRM_ADOPT<Population1, Population2, Solution>::is_feasible_population(const Real tar, const Population & pop) {

		size_t phase = CAST_RP_EPANET->phase();
		size_t interval = CAST_RP_EPANET->interval();
		size_t num = phase*interval;
		Real temp = 0;

		for (size_t i = 0; i < num; ++i) {
			for (int j = 0; j < CAST_RP_EPANET->num_sensor(); ++j) {
				temp += pow(CAST_RP_EPANET->observation_concentration()[j][i], 2);
			}
		}

		Real benchmark = tar * sqrt(temp / (CAST_RP_EPANET->num_sensor()*num));

		size_t count_feasible = 0, count_infeasible = 0;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].objective()[0] <= benchmark) ++count_feasible;
			else ++count_infeasible;
		}

		return count_feasible >= count_infeasible;
	}

	template <typename Population1, typename Population2, typename Solution>
	template <typename Population>
	void LRM_ADOPT<Population1, Population2, Solution>::local_search(Population & pop) {  // 0.2 0.5 0.8
		std::vector<Real> beta = { 0.8f, 0.5f, 0.2f };
		std::vector<std::pair<int, Real>> order_pop;
		int count = -1;
		for (size_t i = 0; i < pop.size(); ++i)
			order_pop.push_back(std::make_pair(++count, pop[i].objective()[0]));
		std::sort(order_pop.begin(), order_pop.end(), CmpByValue());
		int size = order_pop.size();
		int best_idx = order_pop[0].first;
		int sec_best_idx = order_pop[1].first;
		std::vector<int> worst;
		for (size_t i = 0; i < 3; ++i) {
			worst.push_back(order_pop[size - 1 - i].first);
		}

		for (size_t i = 0; i < 3; ++i) {
			pop[worst[i]].variable().index() = pop[best_idx].variable().index();
			pop[worst[i]].variable().start_time() = pop[best_idx].variable().start_time() + beta[i] * (pop[best_idx].variable().start_time() - pop[sec_best_idx].variable().start_time());
			if (pop[worst[i]].variable().start_time() > CAST_RP_EPANET->m_max_start_time)
				pop[worst[i]].variable().start_time() = CAST_RP_EPANET->m_max_start_time;
			else if (pop[worst[i]].variable().start_time() < CAST_RP_EPANET->m_min_start_time)
				pop[worst[i]].variable().start_time() = CAST_RP_EPANET->m_min_start_time;
			for (size_t j = 0; j < pop[worst[i]].variable().multiplier().size(); ++j) {
				pop[worst[i]].variable().multiplier()[j] = pop[best_idx].variable().multiplier()[j] + beta[i] * (pop[best_idx].variable().multiplier()[j] - pop[sec_best_idx].variable().multiplier()[j]);
				if (pop[worst[i]].variable().multiplier()[j] > CAST_RP_EPANET->m_max_multiplier)
					pop[worst[i]].variable().multiplier()[j] = CAST_RP_EPANET->m_max_multiplier;
				else if (pop[worst[i]].variable().multiplier()[j] < CAST_RP_EPANET->m_min_multiplier)
					pop[worst[i]].variable().multiplier()[j] = CAST_RP_EPANET->m_min_multiplier;
			}
		}

	}

}
#endif // OFEC_LRM_ADOPT_LS_H

