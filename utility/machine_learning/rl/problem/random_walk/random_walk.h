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
*  Random_walk problem in RL is a type of classical problem for RL Alg. to evalute the state values.
*
*********************************************************************************/
#ifndef  RL_RANDOM_WALK_H
#define  RL_RANDOM_WALK_H

#include"../problem_base.h"
#include"../../property/property.h"
#include<vector>
#include<string>
#include<map>
#include<algorithm>
#include<iostream>
#include<numeric>

namespace rl {

	class random_walk : public problem_base<RANDOM_WALK_STATE_TYPE, RANDOM_WALK_ACTION_TYPE> {
	public:

		random_walk(int states_num = 11);
		~random_walk();
		void initialize()override;
		RANDOM_WALK_STATE_TYPE & start() { return m_start; }
		const RANDOM_WALK_STATE_TYPE & start() const { return m_start; }
		std::vector<RANDOM_WALK_STATE_TYPE> & terminal() { return m_terminal; }
		const std::vector<RANDOM_WALK_STATE_TYPE> & terminal()const { return m_terminal; }
		std::pair<RANDOM_WALK_STATE_TYPE, double> next_state(const RANDOM_WALK_STATE_TYPE & state, const RANDOM_WALK_ACTION_TYPE & action)override;


	private:
		int m_states_num;
		RANDOM_WALK_STATE_TYPE m_start;
		RANDOM_WALK_STATE_TYPE m_left_end;
		RANDOM_WALK_STATE_TYPE m_right_end;
		std::vector<RANDOM_WALK_STATE_TYPE> m_terminal;
	};


}


#endif // ! RANDOM_WALK_H
