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
*  Td_state_selection template can serve as TD type Algs. with selection state selection.
*  Especially used in situation where actions incur specific   subsequent states.
*
*********************************************************************************/
#ifndef RL_TD_STATE_SELECTION_H
#define RL_TD_STATE_SELECTION_H

#include"./td_update.h"


namespace rl {
	template<typename state_type, typename action_type>
	class td_state_selection :public td_update<state_type, action_type> {

	public:
		td_state_selection() :td_update<state_type, action_type>() {}
		td_state_selection(OFEC::param_map & v) :td_update<state_type, action_type>(v) {}
		virtual ~td_state_selection() {}

		//all states with best state values in limited states;
		std::vector<state_type> best_state(const std::vector<state_type> & states)const {
			std::vector<state_type> beststates;
			float max_statevalue = std::numeric_limits<float>::lowest();
			for (const auto & ele : states) {
				if (this->m_v_value.at(ele).m_state_value > max_statevalue) {
					max_statevalue = this->m_v_value.at(ele).m_state_value;
					beststates.clear();
					beststates.emplace_back(ele);
				}
				else if (this->m_v_value.at(ele).m_state_value == max_statevalue) {
					beststates.emplace_back(ele);
				}
				
			}
			return beststates;
		}

		// randomly choose from the limited states;
		state_type random_state(const std::vector<state_type> & states)const {
			std::uniform_int_distribution<int> uni_int(0, (int)states.size() - 1);
			return states[uni_int(this->m_random_engine)];
		}

		//randomly choose the states with best value in limited states;
		state_type random_best_state(const std::vector<state_type> & states)const {
			std::vector<state_type> all_best = best_state(states);
			return random_state(all_best);
		}

		//In limited states, eplison possibility to choose randomly; 
		//1-epsilon possibility to choose the state with best values;
		state_type epsilon_state(const std::vector<state_type> & states, float epsilon)const {
			std::uniform_real_distribution<float> d_rand;

			if (d_rand(this->m_random_engine) < epsilon) {
				return random_state(states);
				//uniform_int_distribution<int> uni_int(0, actions.size() - 1);
				//return actions[uni_int(m_random_engine)];
			}
			else {
				return random_best_state(states);
			}

		}

	protected:


	private:

	};
}

#endif
