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
// updated Apr 13, 2019 by Li Zhou

#ifndef OFEC_ESs_POPULATION_H
#define OFEC_ESs_POPULATION_H

#include "ESs_Solution.h"
#include "../../../../../../core/algorithm/population.h"


namespace ofec {

	class ESs_population : public population<ESs_Solution>
	{
		using PAIR = std::pair<variable_epanet, Real>;
		struct CmpByValue {
			bool operator()(const PAIR& lhs, const PAIR& rhs) {
				return lhs.second < rhs.second;
			}
		};
	protected:
		std::vector<std::pair<Real, size_t>> m_node_data_obj;
		std::vector<Real> m_probability;
		Real m_delta_user;
		std::vector<std::pair<variable_epanet, Real>> m_whole_pop;
		int m_size_offspring;
		Real m_upper_quantile;
		//int m_iteration = 0;

	public:
		ESs_population(size_t no, size_t dim, size_t size_subpop, Real alpha, Real beta);
		int evolve();
	protected:
		Real cal_standard_deviation(const std::vector<float>& vec1, const std::vector<float>& vec2);
		int produce_offspring(int size_new, Real detla_user);
		int select(Real upper_quantile);

		void update_probability();

		/*template<typename Fun, typename ... Args>
		void sort(Fun fun, Args&& ... args) = delete;*/
	};

}

#endif // !OFEC_ESs_POPULATION_H

