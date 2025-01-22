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

#ifndef OFEC_EVALUATION_DYNAMIC_ENVIRONMENT_H
#define OFEC_EVALUATION_DYNAMIC_ENVIRONMENT_H

#include "../../../../core/environment/environment.h"

namespace ofec {
	class EvaluationDynamicEnvironment : virtual public Environment {
		OFEC_CONCRETE_INSTANCE(EvaluationDynamicEnvironment)
	protected:
		int m_frequency = 0;
		bool m_flag_objective_memory_change = false;
		bool m_flag_variable_memory_change = false;

		void addInputParameters();
		int evaluate_(SolutionBase &sol) override;
		int updateEvaluationTag(SolutionBase &s, Algorithm *alg);
		void handleEvaluateTag(int tag);
	};
}

#endif //!OFEC_EVALUATION_DYNAMIC_ENVIRONMENT_H
