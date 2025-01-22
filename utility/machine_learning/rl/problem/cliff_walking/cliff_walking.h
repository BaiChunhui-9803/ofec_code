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
*  Cliff_walking problem in RL is drived from grid problem, existing cliff reward in the environment.
*
*********************************************************************************/
#ifndef RL_CLIFF_WALKING_H
#define RL_CLIFF_WALKING_H

#include"../grid/grid.h"

namespace rl {
	class cliff_walking :public grid {
	public:
		cliff_walking();
		~cliff_walking();
		//interface to member;
		const GRID_STATE_TYPE & goal() const { return m_goal; }
		GRID_STATE_TYPE & goal() { return m_goal; }
		const GRID_STATE_TYPE & start() const { return m_start; }
		GRID_STATE_TYPE & start() { return m_start; }
		void initialize()override;
		pair<GRID_STATE_TYPE, double>  next_state(const GRID_STATE_TYPE  & state, const GRID_ACTION_TYPE & action)override;

	private:
		//the start&goal point of the grid;
		GRID_STATE_TYPE m_start;
		GRID_STATE_TYPE m_goal;
		//cliff position
		set<GRID_STATE_TYPE> m_cliff_set;
	};

}
#endif // !CLIFF_WALKING
