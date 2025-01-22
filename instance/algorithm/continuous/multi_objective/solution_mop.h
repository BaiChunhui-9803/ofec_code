/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// learning from historical data in the search space for MOPs
// 
//  Created: 28 Feb 2022 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// Last modified: 

#ifndef OFEC_SOLUTION_MOP_H
#define OFEC_SOLUTION_MOP_H


#include "../../../../core/problem/solution.h"

namespace ofec {
	class SolutionMOP :public Solution<> {
	private:
		size_t m_front_type = 0;//0:front sol in subspace
	public:
		size_t getSolType() { return m_front_type; }
		void setSolType(size_t v) { m_front_type = v; }
	};
}

#endif