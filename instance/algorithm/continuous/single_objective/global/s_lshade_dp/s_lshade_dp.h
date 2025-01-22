/********* Begin Register Information **********
{
	"name": "S_LSHADE_DP",
	"identifier": "S_LSHADE_DP",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------/
*
*********************************************************************************/

#ifndef OFEC_S_LSHADE_DP_H
#define OFEC_S_LSHADE_DP_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class S_LSHADE_DP : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(S_LSHADE_DP)
	protected:
		using Individual = ofec::Continuous::SolutionType;

		
		void addInputParameters();
		void initialize_(Environment* env) override;
		void run_(Environment* env) override;
		

		void modifySolutionWithParentMedium(Individual& child, const Individual& parent, Environment* env);
		void reducePopulationWithSort(std::vector<Individual>& pop, std::vector<int>& stagnation);


		void operateCurrentToPBest1BinWithArchive(const std::vector<Individual>& pop, Individual& child, int& target, int& p_best_individual, double& scaling_factor, double& cross_rate, const std::vector<Individual>& archive, int& arc_ind_count, Environment* env);
		void operateTarget1BinWithArchive(const std::vector<Individual>& pop, Individual& child, int& target, double& scaling_factor, double& cross_rate, const std::vector<Individual>& archive, int& arc_ind_count, Environment* env);

		void updateBsf(const std::vector<Individual>& pop, Individual& bsf);
		void updateDatum(const std::vector<Individual>& pop,Environment* env);
	protected:
		//L-SHADE parameters
		size_t m_pop_size = 100;
		size_t m_memory_size = 6;
		double m_arc_rate = 2.6;
		double m_p_best_rate = 0.11;
		int reduction_ind_num = 0;

		int m_arc_size = 0;

		int NG = 20;
		int M = 2;
		double gamma = 0.3;

		std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)> m_eval_fun;

	};
}

#endif