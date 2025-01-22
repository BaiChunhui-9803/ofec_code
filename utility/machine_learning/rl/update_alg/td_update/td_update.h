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
*  TD_update template can serve as all TD type Algs. with state_type and aciton type fixed.
*
*********************************************************************************/
#ifndef RL_TD_UPDATE_H
#define RL_TD_UPDATE_H

#include"../updata_base.h"

namespace rl {
	template<typename state_type, typename action_type>
	class td_update :public update_base<state_type, action_type> {

	public:
		td_update() :update_base<state_type, action_type>(), m_discount_factor(0.8) {}
		td_update(OFEC::param_map & v) :update_base<state_type, action_type>(v), m_discount_factor(v.at("discountFactor")) {}
		virtual ~td_update() {}

		virtual void initialize(const std::unordered_map<state_type, std::vector<action_type> > & state_actions, \
			float initial_v = 0.f, float initial_v_error = 0.f, float initial_q = 0.f, float initial_q_error = 0.f, \
			bool if_set_actions = true) {
			update_base<state_type, action_type>::initialize(state_actions, initial_v, initial_v_error, initial_q, initial_q_error, if_set_actions);
		}

		// on_policy TD(0) control alg. update (Sarsa) with overloading the function: update_q;
		virtual void update_q(const state_type & state, const action_type & action, const state_type & next_state,
			const action_type & next_action, float r) {
			this->update_q_(state, action, next_state, next_action, r);

		}

		// off_policy TD(0) control alg. update (Q-learning) with overloading the function: update_q;
		virtual void update_q(const state_type & state, const action_type & action, const state_type & next_state, float r) {
			this->update_q_(state, action, next_state, r);
		}

		// TD(0) predication alg. for state_value;
		virtual void update_v(const state_type & state, const state_type & next_state, float r) {
			//this->m_v_value[state] +=0.01*(r + m_discount_factor*this->m_v_value[next_state] - this->m_v_value[state]);
			this->update_v_(state, next_state, r);
		}

		//interface of the dicount factor;
		float discount_factor()const { return m_discount_factor; }
		void set_discount_factor(float discount_factor) { m_discount_factor = discount_factor; }

	protected:
		float m_discount_factor;
	private:

		void update_q_(const state_type & state, const action_type & action, const state_type & next_state, float r) {
			action_type max_action = update_base<state_type, action_type>::random_best_action(next_state);
			//if (!std::isnormal(this->m_q_value[std::make_pair(state, action)])) {
			//	std::cout << "bug triggered" << std::endl;
			//}
			//auto value_before_bug = this->m_q_value[make_pair(state, action)];
			action_value_type & modified_data = this->m_q_value[std::make_pair(state, action)];
			float estimation_error = r + m_discount_factor * this->m_q_value[std::make_pair(next_state, max_action)].m_action_value - modified_data.m_action_value;

			modified_data.m_error = (1 - this->m_learning_ratio)* modified_data.m_error + this->m_learning_ratio*std::fabs(estimation_error);
			modified_data.m_action_value += this->m_learning_ratio*estimation_error;
			//	(r + m_discount_factor * this->m_q_value[make_pair(next_state, max_action)] - this->m_q_value[make_pair(state, action)]);
			//if (!std::isnormal(this->m_q_value[std::make_pair(state, action)])) {
			//	std::cout << "bug triggered" <<"  and the value before the bug is "<<value_before_bug<< std::endl;
			//}
		}
		void update_v_(const state_type & state, const state_type & next_state, float r) {
			//this->m_v_value[state] +=0.01*(r + m_discount_factor*this->m_v_value[next_state] - this->m_v_value[state]);
			state_value_type & modified_data = this->m_v_value[state];

			float estimation_error = r + m_discount_factor * this->m_v_value[next_state].m_state_value - modified_data.m_state_value;
			modified_data.m_error = (1 - this->m_learning_ratio)*modified_data.m_error + this->m_learning_ratio*std::fabs(estimation_error);
			modified_data.m_state_value += this->m_learning_ratio* estimation_error;
			//this->m_v_value[state] += this->m_learning_ratio*(r + m_discount_factor * this->m_v_value[next_state] - this->m_v_value[state]);
		}
		void update_q_(const state_type & state, const action_type & action, const state_type & next_state,
			const action_type & next_aciton, float r) {
			action_value_type & modified_data = this->m_q_value[std::make_pair(state, action)];
			float estimation_error = r + m_discount_factor * this->m_q_value[std::make_pair(next_state, next_aciton)].m_action_value - modified_data.m_action_value;
			modified_data.m_error = (1 - this->m_learning_ratio)*modified_data.m_error + this->m_learning_ratio*std::fabs(estimation_error);
			//this->m_q_value_error[std::make_pair(state, action)] = (1 - this->m_learning_ratio)*this->m_q_value_error[std::make_pair(state, action)] + this->m_learning_ratio*std::fabs(estimation_error);
			//class template can resolve the elements of its superclass by this ptr;
			modified_data.m_action_value += this->m_learning_ratio*estimation_error;
			//	(r + m_discount_factor * this->m_q_value[make_pair(next_state, next_aciton)] - this->m_q_value[make_pair(state, action)]);
		}
	};
}

#endif // !Q_LEARNING_H