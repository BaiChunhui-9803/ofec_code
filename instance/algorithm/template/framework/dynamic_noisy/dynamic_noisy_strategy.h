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
* class Algorithm is an abstract for all algorithms.
*
*********************************************************************************/
#ifndef OFEC_DYNAMIC_NOISY_STRATEGY_H
#define OFEC_DYNAMIC_NOISY_STRATEGY_H
#include"../../../../../core/definition.h"
#include<vector>

namespace ofec {

	
	class DynamicNoisyStartegy  {
	protected:
		Real m_uncertianty_thresHold = 0.005;
	public:
		DynamicNoisyStartegy() = default;
		virtual ~DynamicNoisyStartegy() = default;
		bool judgeChange(Problem *pro, const std::vector<Real>& originObjs, const std::vector<Real>& curObjs) const;

	};
}

#endif