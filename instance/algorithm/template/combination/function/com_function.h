
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
* some general methods used in ofec are defined in this file, including distance
* metrics, solution domination relationship,
*
*********************************************************************************/

#ifndef OFEC_COM_FINCTIONAL_H
#define OFEC_COM_FINCTIONAL_H

#include "../../../../../core/definition.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {

	//template<typename TVecIndi>
	//inline void updateFitness(Problem *pro,const TVecIndi& pop)
	//{
	//	
	//	if (pro->hasTag(ProblemTag::MOP)) {
	//		std::vector<std::vector<Real>*> data(pop.size());
	//		for (size_t i = 0; i < pop.size(); ++i)
	//			data[i] = &(m_individuals[i].objective());
	//		std::vector<int> rank;
	//		//nd_sort::fast_sort<Real>(data, rank, global::ms_global->m_problem->optimizeMode());
	//		if (m_individuals.size() > 1e3) {
	//			nd_sort::filter_sort_p<Real>(data, rank, pro->optimizeMode());
	//		}
	//		else {
	//			nd_sort::filter_sort<Real>(data, rank, pro->optimizeMode());
	//		}
	//		sort(pro);
	//		for (auto& it : m_individuals) {
	//			it->setFitness(it->fitness());
	//		}
	//	}
	//	else {
	//		if (pro->optimizeMode()[0] == OptimizeMode::Maximize) {
	//			for (auto& it : m_individuals) {
	//				it->setFitness(1.0 / (it->objective(0) + 1e-5));
	//			}
	//		}
	//		else {
	//			for (auto& it : m_individuals) {
	//				it->setFitness(it->objective(0));
	//			}
	//		}
	//	}
	//}

	//template<typename TInd>
	//inline void updateFitness(Problem *pro, TInd& sol) const
	//{
	//	if (pro->hasTag(ProblemTag::MOP)) {
	//		sol.setFitness(sol.fitness());
	//	}
	//	else {
	//		if (pro->optimizeMode()[0] == OptimizeMode::Maximize) {
	//			sol.setFitness(1.0 / (1e-5 + sol.objective(0)));
	//		}
	//		else {
	//			sol.setFitness(sol.objective(0));
	//		}
	//	}
	//}



}

#endif