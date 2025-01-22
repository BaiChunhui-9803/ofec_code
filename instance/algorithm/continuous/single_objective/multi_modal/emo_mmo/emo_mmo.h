/********* Begin Register Information **********
{
	"name": "EMO-MMO",
	"identifier": "EMO_MMO",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
@article{cheng2017evolutionary,
  volume = {22},
  number = {5},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {692--706},
  year = {2017},
  author = {Ran Cheng and Miqing Li and Ke Li and Xin Yao},
  title = {Evolutionary multiobjective optimization-based multimodal optimization: Fitness landscape approximation and peak detection}
}
*-------------------------------------------------------------------------------
* Implemented by Junchen Wang (wangjunchen.chris@gmail.com) at 2020/9/26
*********************************************************************************/

#ifndef OFEC_EMO_MMO_H
#define OFEC_EMO_MMO_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../template/classic/genetic_algorithm/sbx_pop.h"
#include "../../../../template/classic/differential_evolution/population.h"
#include <list>
#include <set>
#include <unordered_set>

namespace ofec {
	class transferedMOP {
	public:
		transferedMOP(size_t t_max, size_t size_var, size_t size_pop);
		void updateNormX(const Population<Solution<>> &offspring);
		void updateDelta(int t);
		void evaluate2ndObj(Population<Solution<>> &offspring);
	private:
		void updateManhDis();
		int minManhDis(size_t i);
		size_t m_t_max, m_N, m_size_var;
		Real m_delta;
		std::vector<std::vector<int>> m_norm_x, m_manh_dis;
	};

	class EMO_MMO : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(EMO_MMO)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	private:
		void MOFLA(Environment *env);
		void binaryCuttingAPD(Environment *env);
		void localSearch(Environment *env);
		void landscapeCutting(std::unordered_set<size_t> &D_c, Real eta, Environment *env);
		void APD(std::vector<std::vector<Real>> &manh_dis_mat, std::unordered_set<size_t> D_c, Environment *env);
		void updateManhDisMat(std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c);
		void normcManhDisMat(std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c, Environment *env);
		Real minManhDis(const std::vector<std::vector<Real>> &manh_dis_mat, const std::unordered_set<size_t> &D_c, size_t i);
		size_t m_size_pop, m_t_max;
		Real m_ratio_MOFLA_FEs;                         // percentage of FEs allocated to fitness landscape approximation
		Real m_eta;                                     // initial cutting ratio
		std::unique_ptr<PopSBX<>> m_MOFLA_pop;                          // population for MOFLA
		MultiPopulation<PopulationDE<>> m_LS_pops;   // multipopulation for local search
		std::vector<std::unique_ptr<Solution<>>> m_D;   // external archive for approximating fitness landscape
		std::set<std::vector<Real>> m_D_x;              // decision variables of the external archive for uniqueness
		std::vector<std::list<size_t>> m_P;             // detected peak set
	};
}

#endif // !OFEC_EMO_MMO_H

