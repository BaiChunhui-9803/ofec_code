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
* class base defines the virtual base for model, problem and update alg. in reinforcement learning.
*
*********************************************************************************/
#ifndef RL_BASE_H
#define RL_BASE_H

namespace rl {
	

	class base {

	public:

		base() = default;
		virtual ~base() {};
		base(const base&) = default;
		base(base&&) = default;
		base& operator=(const base&) = default;
		base& operator=(base&&) = default;
	protected:

	private:

	};
}



#endif // !BASE_H