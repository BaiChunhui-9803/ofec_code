/******************************************************************************
* Project: Framework for Reinforcement Learning (Also the part of utilities in OFEC)
*******************************************************************************
* Author: Xia Hai
* Email: strawberry9583@gmail.com
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
* class template simle_model serves for RL experience storage with table as an example of the usage of model_base.
*
*********************************************************************************/
#ifndef RL_SIMPLE_MODEL_H
#define RL_SIMPLE_MODEL_H

#include<algorithm>
#include"model_base.h"
#include"../property/property.h"

namespace rl {
	template<typename state_type, typename action_type>
	class simple_model :public model_base<state_type, action_type> {
	public:
		simple_model() :model_base<state_type, action_type>() {};
		~simple_model() {};

		void feed(const state_type & state, const action_type & action, const state_type & next_state, double reward) {
			// S_A pair hasn't occurred.
			if (std::find(m_int2sa.begin(), m_int2sa.end(), std::make_pair(state, action)) == m_int2sa.end()) {
				m_int2sa.push_back(std::make_pair(state, action));
			}
			//add the SASR sequence into model;
			m_model_data[std::make_pair(state, action)].push_back(std::make_pair(next_state, reward));
		};

		std::tuple<state_type, action_type, state_type, double> sample() {
			std::uniform_int_distribution<int> rand(0, (int)m_int2sa.size() - 1);
			int s_a_idx = rand(this->m_random_engine);
			std::uniform_int_distribution<int> rand1(0, (int)m_model_data[m_int2sa[s_a_idx]].size() - 1);
			int s_r_idx = rand1(this->m_random_engine);
			return make_tuple(m_int2sa[s_a_idx].first, m_int2sa[s_a_idx].second,
				m_model_data[m_int2sa[s_a_idx]][s_r_idx].first,
				m_model_data[m_int2sa[s_a_idx]][s_r_idx].second);
		};
	private:
		// record the S_A pair happend; though also stored in m_model_data, but will speed up sample;
		// because the map can not be randomly accessed.
		std::vector<std::pair<state_type, action_type>> m_int2sa;
		//record the SASR sequences;
		std::map<std::pair<state_type, action_type>, std::vector<std::pair<state_type, double>>> m_model_data;
	};
}

#endif // !SIMPLE_MODEL_H
