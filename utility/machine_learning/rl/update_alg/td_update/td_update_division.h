///******************************************************************************
//* Project: Framework for Reinforcement Learning (Also the part of utilities in OFEC)
//*******************************************************************************
//* Author: Xia Hai
//* Email: strawberry9583@gmail.com
//* Language: C++
//*-------------------------------------------------------------------------------
//*  This file is part of OFEC. This library is free software;
//*  you can redistribute it and/or modify it under the terms of the
//*  GNU General Public License as published by the Free Software
//*  Foundation; either version 2, or (at your option) any later version.
//*
//*  see https://github.com/Changhe160/OFEC for more information
//*
//*-------------------------------------------------------------------------------
//*  td_update_division class serves as the subsidiary of the reinforcement learning based evolutionary computation 
//*  with solution space division.
//*
//*********************************************************************************/
//#ifndef RL_TD_UPDATE_DIVISION_H
//#define RL_TD_UPDATE_DIVISION_H
//
//#include"./td_update.h"
//#include<list>
//
//namespace rl {
//	class td_update_division:public td_update<int,int> {
//
//
//	public:
//		td_update_division(OFEC::param_map & v);
//		~td_update_division ();
//		const std::unordered_map<std::pair<int, int>, int,rl::state_action_pair_hash> & visit_count() const;
//		std::unordered_map<std::pair<int, int>, int,rl::state_action_pair_hash> & visit_count();
//		void initialize(const std::unordered_map<int, std::vector<int> > & state_actions, \
//			float initial_v = 0.f, float initial_v_error = 0.f, float initial_q = 0.f, float initial_q_error = 0.f)override;
//		//add new state_action value;
//		void add_action_value(std::unordered_map<std::pair<int,int>, rl::action_value_type,rl::state_action_pair_hash> & added_value);
//		//splite the state into two halves with q_value and v_value added too;
//		void splite_state(int state, const std::vector<int> & neighbor_link);
//		//override the function update_q;
//		void update_q(const int & state, const int & action, const int & next_state, float r)override;
//	private:
//		std::unordered_map<std::pair<int,int>,int,rl::state_action_pair_hash> m_visit_count;//count the sampling times of each possible state_action pair;
//	};
//
//
//}
//
//#endif // !RL_TD_UPDATE_DIVISION_H
