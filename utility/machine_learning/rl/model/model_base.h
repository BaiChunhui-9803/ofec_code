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
* class template model_base defines the virtual base for model to store the RL experience.
*
*********************************************************************************/
#ifndef RL_MODEL_BASE_H
#define RL_MODEL_BASE_H

#include"base.h"
#include"../property/property.h"
#include<map>
#include<vector>
#include<random>

namespace rl {
	template<typename state_type, typename action_type>
	class model_base :public base {

	public:

		model_base() { m_random_engine.seed(0); };
		virtual ~model_base() {};
		virtual void feed(const state_type & state, const action_type & action, const state_type & next_state, double reward) {}
		virtual std::tuple<state_type, action_type, state_type, double> sample() {
			return std::make_tuple(state_type(), action_type(), state_type(), double());
		}

	protected:

		//mutable (for the const function access safety) random engine for data storage model;
		mutable std::default_random_engine m_random_engine;

	};

}

#endif // !MODEL_BASE_H