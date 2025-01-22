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

#ifndef OFEC_SELECTION_METHODS_H
#define OFEC_SELECTION_METHODS_H

#include <cmath>
#include <vector>
#include "../../core/definition.h"
#include "../../utility/nondominated_sorting/fast_sort.h"
#include "../../utility/nondominated_sorting/filter_sort.h"
#include "../../utility/functional.h"

namespace ofec {
	//NSGA-II environmental selection
	template<typename TPopulation>
	void popNondominatedSorting(TPopulation& pop, std::vector<OptimizeMode> m_optimize_mode) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i)
			objs.emplace_back(&pop[i].objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
		for (size_t i = 0; i < pop.size(); ++i)
			pop[i].setFitness(rank[i]);
	}

	template<typename T>
	std::vector<int> vectorNondominatedSorting(std::vector<std::vector<T>>& data, std::vector<OptimizeMode> m_optimize_mode) {
		std::vector<std::vector<T>*> objs;
		for (size_t i = 0; i < data.size(); ++i)
			objs.emplace_back(&data[i]);
		std::vector<int> rank;
		nd_sort::fastSort<T>(objs, rank, m_optimize_mode);
		return rank;
	}

	template<typename TPopulation, typename TPopCombined>
	void envirSelectionByNSGAII(TPopulation& pop, TPopCombined& pop_combined, std::vector<OptimizeMode> m_optimize_mode) {
		popNondominatedSorting(pop_combined, m_optimize_mode);
		int pops = 0;  //indicate parent population size be 0
		int size = pop_combined.size();
		int rank = 0;
		size_t m_num_obj = m_optimize_mode.size();
		while (true) {
			int count = 0;
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					count++;
			int size2 = pops + count;
			if (size2 > pop.size()) {
				break;
			}
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
				{
					pop[pops] = pop_combined[i];
					++pops;
				}
			rank++;
			if (pops >= pop.size()) break;
		}
		if (pops < pop.size()) {
			std::vector<int> list;
			// save the Solutions in the overflowed front
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					list.push_back(i);
			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < m_num_obj; j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop_combined[list[i]].objective()[j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = pop.size() - pops;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				pop[pops] = pop_combined[list[idx[i]]];
				++pops;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
	}

	//select from the first layer
	template<typename T>
	void crowSelectInFirstLayer(std::vector<std::vector<T>>& vector1, std::vector<std::vector<T>>& vector2) {
		//vector2 is the first layer vector
		int pops = 0;  //indicate parent population size be 0
		int size = vector2.size();
		int rank = 0;
		size_t m_num_obj = vector1[0].size();
		while (true) {
			int count = vector2.size();
			int size2 = pops + count;
			if (size2 > vector1.size()) {
				break;
			}
		}

		if (pops < vector1.size()) {
			std::vector<int> list;
			// save the Solutions in the overflowed front
			for (size_t i = 0; i < size; i++) {
				list.push_back(i);
			}

			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < m_num_obj; j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = vector2[list[i]][j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = vector1.size() - pops;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				vector1[pops] = vector2[list[idx[i]]];
				++pops;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
	}
}

#endif // !OFEC_SELECTION_METHODS_H