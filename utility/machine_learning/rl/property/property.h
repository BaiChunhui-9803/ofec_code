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
*  Property defines the popular state_type and action_type in RL.
*
*********************************************************************************/
#ifndef RL_PROPERTY_H
#define RL_PROPERTY_H

#include<tuple>

namespace rl {


	enum GRID_ACTION_TYPE {
		GRID_ACTION_UP, GRID_ACTION_DOWN, GRID_ACTION_LEFT, GRID_ACTION_RIGHT
	};
	enum RANDOM_WALK_ACTION_TYPE {
		RANDOM_WALK_ACTION_LEFT, RANDOM_WALK_ACTION_RIGHT
	};


	typedef std::pair<int, int> GRID_STATE_TYPE;
	typedef int RANDOM_WALK_STATE_TYPE;

	typedef std::tuple<GRID_STATE_TYPE, GRID_ACTION_TYPE, GRID_STATE_TYPE, double> GRID_SASR_QUADRUPLE;
	typedef std::tuple<RANDOM_WALK_STATE_TYPE, RANDOM_WALK_ACTION_TYPE, RANDOM_WALK_STATE_TYPE, double> RANDOM_WALK_SASR_QUADRUPLE;

	
}


#endif // !DEFINITION_H
