/********* Begin Register Information **********
{
	"name": "KDT_MOEA_DE",
	"identifier": "KDT_MOEA_DE",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Xiaofang Wu
* Email: changhe.lw@google.com Or wuxiaofang@cug.edu.cn
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 11 July 2019
// Last modified: 12 Sep 2019 by Xiaofang Wu (email:email:wuxiaofang@cug.edu.cn)

#ifndef OFEC_KDT_MOEA_DE_H
#define OFEC_KDT_MOEA_DE_H

#include "../kdt_moea.h"
#include "../../moea_de_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../predict_PF.h"
#include "../../../../template/selection/multi_objective/moead.h"

namespace ofec {
	class IndDE;
	class KDT_MOEA_DE_pop:public PopMODE<IndDE>{
	public:
		KDT_MOEA_DE_pop(size_t size_pop, Problem *pro,const ParameterMap& v);
		void initialize_(Problem *pro, Random *rnd);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
	protected:
		int evolve_mo(Problem *pro, Algorithm *alg, Random *rnd);
	public:
		/***********************MOEA/D-start*************************/
		int evolve_MOEAD(PopMODE<>& pop_MOEAD,Problem *pro, Random *rnd);
		int evolve_Mo_MOEAD(PopMODE<>& pop_MOEAD,Problem *pro, Random *rnd);
		void initialize_MOEAD(std::vector<std::unique_ptr<IndDE>>& parent,Problem *pro);
		void update_reference_MOEAD(IndDE& sol,Problem *pro);
		void init_uniformweight_MOEAD(int parent_size, Problem *pro);
		void min_max_sort(std::vector<Real>& v, std::vector<int>& idx);
		void init_neighbourhood_MOEAD();
		void update_problem_MOEAD(std::vector<std::unique_ptr<IndDE>>& parent, IndDE& sol, int id, int type,Problem *pro, Random *rnd);
		Real fitnessfunction_MOEAD(std::vector<Real>& obj, int k,Problem *pro);
		//void update_reference(DE::Solution& sol);
		//void update_problem_MOEAD(std::vector<std::unique_ptr<DE::Solution>>& parent, DE::Solution& sol, int id, int type);
		std::vector<Real> mv_ideal_point;	//the best value in every dimension
		std::vector<IndDE> mv_sol_arr;	//corresponding Solution of the best value in every dimension
		std::vector<Vector> mv_namda;
		std::vector<std::vector<int> > mvv_neigh;
		int m_limit = 2;
		int m_niche = 20;		//number of neighbours
		int m_unit = 12;
		std::vector<int> m_type;
		enum DecomFun { _TCHE1, _TCHE2, _NBI1, _NBI2, _NBI3 };
		DecomFun m_decom_function = _TCHE1;
		/***********************MOEA/D-end***************************/

		/*EvaluationTag evolve();
		EvaluationTag evolve_mo();*/
		int evolve_mo_new(Problem *pro, Algorithm *alg, Random *rnd);
		//void eval_selection();
		//void eval_selection_();		//select an Solution every cycle
		//void eval_selection__();	//select an Solution every cycle, and consider F1 of nondominated-sorting
		void eval_selection_new(Problem *pro);	//select an Solution every cycle, and consider F1 of nondominated-sorting, and sample reference point, and not use cdi indicator

		int find_non_dominance(std::vector<IndDE>& offspring,Problem *pro);
		void select_F2_Fn_new(std::vector<IndDE>& F2_Fn, int& pops,Problem *pro);
		void select_F2_Fn(std::vector<IndDE>& F2_Fn, int& pops,Problem *pro);
		void select_F1(std::vector<IndDE>& F1, int& pops,Problem *pro);
		int update_region_to_sort_ind(std::vector<IndDE>& offspring, std::vector<std::vector<int>>& region_to_sort_ind);
		void select_an_ind_subspace_center(int& pops, std::vector<IndDE>& F1, std::vector<std::vector<int>>& region_to_sort_ind, int k,Problem *pro);//select an Solution from each subapces, let them near the center of subspace
		void select_M_ind_F1(int& pops, std::vector<IndDE>& F1,Problem *pro);//select M Solution from F1, let them near the objective axis

		void select_an_ind_subspace_diversity(int& pops, std::vector<IndDE>& F2_Fn, std::vector<std::vector<int>>& region_to_sort_ind, int k,Problem *pro);//select an Solution from non-dominated subapces, let them keep diversity with selected population


		void eval_selection_dd_dominance(Problem *pro);
		void eval_selection_1_k_dominance(Problem *pro);
		void eval_selection_ranking_dominance(Problem *pro);
		//int dominance_1_k(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<optimization_mode>& opt_mode);
		//dominationship objective_compare_1_k(const std::vector<Real>& a, const std::vector<Real>& b, const std::vector<optimization_mode> &mode);
		int ranking_dominance_R_sum(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<OptimizeMode>& opt_mode);
		//int ranking_dominance_R_min(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<optimization_mode>& opt_mode);

		void find_non_dominance(Problem *pro);	//find non-dominated Solutions from m_pop, and store them into m_arhive
		void update_index_archive();//update index of inviduals according to their distance
		//void prediction(const DE::MOEA_DE_pop<>& Pop, int pops); //predict N Solutions according to the selected-Solutions
		void predict();	//predict some points according to m_archive by kalman filter and least square
		void add_pop_to_pred(Problem *pro); //add the non-dominated Solutions of pre-generation population into prediction (its left lower point in box)
		void add_subregion_to_pred(const std::vector<std::pair<int, Real>>& region_dds, int type,Problem *pro); //add the non-dominated(type=0) or all(type=1, sub_dds>0) subregions into prediction (its left lower point in box)
		void add_shrink_F1_to_pred(); //add shrinked non-dominated Solutions into prediction
		void add_shrink_all_nondominance_to_pred(const std::vector<IndDE>& all_nondominance); //add shrinked m_all_nondominance into prediction
		void handle_prediction(const std::vector<IndDE>& offspring,Problem *pro);//erase these dominated-prediction-point by offspring

		bool judge_point(const std::vector<Real>& point, const std::vector<std::vector<Real>>& prediction);//judge whether prediction include point

		void nondominatedSorting(std::vector<IndDE>& offspring,Problem *pro);	//NSGAII
		void eval_selection_NSGAII(PopMODE<>& parent, std::vector<IndDE>& offspring,Problem *pro); //NSGAII

		//void find_dd();
		void select_ind_F1(int& pops, std::vector<int>& F1,Problem *pro);
		void select_ind_F1(int& pops, std::vector<IndDE>& F1,Problem *pro);
		void select_ind_subregion(int& pops, std::vector<int>& vec_inds, int num_max,Problem *pro);
		void select_an_ind_max_min(int& pops, std::vector<IndDE>& offspring, std::vector<std::pair<int, std::pair<int, int>>>& sub_inds, std::vector<std::vector<int>>& region_to_sort_ind); //the first store index of Solution, the second store index of subspace
		Real dd(IndDE& ind1, IndDE& ind2,Problem *pro);
		void recordObj(int iter, Algorithm *alg, Problem *pro);//
		//void normalization(const std::vector<DE::Solution> &offspring, std::vector<DE::Solution> &offspring_map, const std::vector<std::pair<Real, Real>> &min_max);	//map m_offfspring to m_offspring_map, i.e., map objective value to [0,1]
		void setAllNondomiSolu() { m_A = m_F1; }
		std::vector<IndDE>& getAllFrontSolu() { return m_F1; }
	protected:
		//int m_not_pop, m_nsga2_pop, m_dd_pop, m_both_pop;
		int m_number_test = 11;	//test
		std::vector<IndDE> m_pop_rank_is_0;	//record Solutions of rank=0		
		std::vector<IndDE> m_pop_first_select; //record Solutions of the first selection, size is count_F1 or M
		std::vector<IndDE> m_pop_second_select;//record Solutions of the second selection, size is 1

		std::vector<IndDE> m_pop_1_k_dominance;	//test
		std::vector<IndDE> m_pop_ranking_dominance;//test

		Real m_k = 0.5;	//parameter for (1-k)-dominance

		//DE::MOEA_DE_pop<> m_pop;
		std::vector<IndDE> m_offspring;	// 2 size of population
		std::vector<IndDE> m_Q;			// size of population, only store Q
		std::unique_ptr<KDT_MOEA<PopMODE<>, IndDE>> m_kdt;

		//std::vector<std::vector<Real>> m_all_nondominance;	//store all non-dominated Solutions in historical Solutions 
		std::vector<IndDE> m_A;	//store all non-dominated Solutions in historical Solutions, A^r record m_F1, A^d record m_F2_Fn 

		std::vector<std::vector<IndDE>> m_archive;	//store non-dominated Solutions every generation
		std::vector<kalmanFilter<>> m_Kalman_Filter;
		leastSquare<> m_least_square;
		std::vector<std::vector<Real>> m_prediction_kalman_filter;
		std::vector<std::vector<Real>> m_prediction_least_square;
		std::vector<std::vector<Real>> m_ref_point;
		std::vector<IndDE> m_prediction;
		int m_prediction_size = -1;
		int m_prediction_handled_size = -1;

		int m_count_F1 = 0;
		std::vector<IndDE> m_F1;
		std::vector<IndDE> m_F2_Fn;
		
	public:
		std::vector<std::vector<std::vector<Real>>> m_observe_angle_subspace_inner;
		bool m_observe = true;
		PopMODE<> m_pop_MOEAD;
	};

	class KDT_MOEA_DE :public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();
#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
		std::unique_ptr<KDT_MOEA_DE_pop> m_pop;
		size_t m_pop_size;
		
	};

}

#endif // !OFEC_KDT_MOEA_DE_H
