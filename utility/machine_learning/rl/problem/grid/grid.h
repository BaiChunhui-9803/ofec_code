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
*  grid type problem in RL is the ususal enviroment with meshes.
*
*********************************************************************************/
#ifndef RL_GRID_H
#define RL_GRID_H

#include"../problem_base.h"
#include"../../property/property.h"
#include<vector>
#include<string>
#include<map>
#include<algorithm>
#include<iostream>
#include<set>

namespace rl {

	using namespace std;

	class grid :public problem_base<GRID_STATE_TYPE, GRID_ACTION_TYPE> {
	public:
		grid(int width = 12, int height = 4);
		virtual ~grid();
		void initialize()override;
		pair<GRID_STATE_TYPE, double>  next_state(const GRID_STATE_TYPE  & state, const GRID_ACTION_TYPE & action)override;

		//void reset_obstacle(const set<pair<int, int>> & obstacle);
		//void print_grid();
		//void print_optimal_policy();
	private:

	protected:
		//the grid size;
		int m_width;
		int m_height;

		//the obstacle of the grid;
		//set<GRID_STATE_TYPE> m_obstacle;

	};
}


#endif // 
