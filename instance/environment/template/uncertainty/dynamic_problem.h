/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
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
*
*********************************************************************************/

// Created: 9 June 2021
// Modified:  2024 04 23 by DIAOYIYA 
// Last modified:

#ifndef OFEC_DYNAMIC_PROBLEM_H
#define OFEC_DYNAMIC_PROBLEM_H

#include "../../../../core/problem/problem.h"

namespace ofec {
	class DynamicProblem : virtual public Problem {
	public:
		virtual void change(Random *rnd) {}
	};
}

#endif