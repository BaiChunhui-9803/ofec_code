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
// Created: 11 Mar 2024


#ifndef OFEC_SPACING_H
#define OFEC_SPACING_H

#include "../functional.h"
#include <vector>
#include <numeric>

namespace ofec {

	template<typename T>
	Real Spacing(std::vector<std::vector<T>>& pop) {
		//std::vector<std::vector<T>> pop_temp = pop;
		Real score = 0.;
		int pop_size = pop.size();
		int obj_size = pop[0].size();
		std::vector<Real> ind_dist;
		Real sum_dist = 0.;
		for (size_t i = 0; i < pop_size; ++i) {
			auto ind = pop[i];
			Real dist = 1.e14;
			for (size_t j = 0; j < pop_size; ++j) {
				if (i != j) {
					auto m_dist = manhattanDistance(ind.begin(), ind.end(), pop[j].begin());
					if (m_dist < dist)
						dist = m_dist;
				}
			}
			sum_dist += dist;
			ind_dist.push_back(dist);
		}
		Real mean_dist = sum_dist / pop_size;
		for (size_t i = 0; i < pop_size; ++i) {
			score += std::pow((mean_dist - ind_dist[i]), 2);
		}
		return std::sqrt(score / (pop_size - 1));
	}
}

#endif
