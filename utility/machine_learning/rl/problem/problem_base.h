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
*  All problems in RL will be derived from the problem base.
*
*********************************************************************************/
#ifndef RL_PROBLEM_BASE_H
#define RL_PROBLEM_BASE_H

#include"../model/base.h"
#include<vector>
#include<unordered_map>
#include<map>

namespace rl {

	template<typename state_type, typename action_type>
	class problem_base : public base {
	public:
		problem_base() = default;
		problem_base(int state_num) : m_state_num(state_num) {}
		int state_num()const { return m_state_num; }
		virtual ~problem_base() {}
		virtual void initialize() {}
		virtual void set_actions4states() {}
		void set_states(std::vector<state_type> & states) { m_states = states; }
		std::vector<state_type> & states() { return m_states; }
		const std::vector<state_type> & states()const { return m_states; }
		void set_actions(std::vector<action_type> & actions) { m_actions = actions; }
		std::vector<action_type> & actions() { return m_actions; }
		const std::vector<action_type> & actions()const { return m_actions; }
		void set_states2actions(std::unordered_map<state_type, std::vector<action_type> > & states2actions) {
			m_states2actions = states2actions;
		}
		std::unordered_map<state_type, std::vector<action_type> > & states2actions() { return m_states2actions; }
		const std::unordered_map<state_type, std::vector<action_type> > & states2actions()const { return m_states2actions; }
		void set_states2reward(std::map<state_type, double> & states2rewards) {
			m_state2reward = states2rewards;
		}
		std::unordered_map<state_type, double> & states2rewards() { return states2rewards; }
		const std::unordered_map<state_type, double> & states2rewards()const { return states2rewards; }

		virtual std::pair<state_type, double> next_state(const state_type & state, const action_type & action) {
			return std::make_pair(state_type(), 0.);
		}

	protected:
		//state num;
		int m_state_num;
		//states set;
		std::vector<state_type> m_states;
		//overall action set;
		std::vector<action_type> m_actions;
		//action set conresponding to the state;
		std::unordered_map<state_type, std::vector<action_type>> m_states2actions;
		//returns conresponding to the the state;
		std::unordered_map<state_type, double> m_state2reward;
	};

}



#endif