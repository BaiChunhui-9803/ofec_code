/********* Begin Register Information **********
{
	"name": "MMOP_CEC2015_F07",
	"identifier": "CEC2015_MMOP_F07",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Problem Definitions and Evaluation Criteria for the CEC 2015
*  Competition on Single Objective Multi-Niche Optimization.
*******************************************************************************************/


#ifndef OFEC_F7_SR_EXPANDED_SIX_HUMP_CAMEL_BACK_H
#define OFEC_F7_SR_EXPANDED_SIX_HUMP_CAMEL_BACK_H

#include "function.h"

namespace ofec {
	namespace cec2015 {
		class F7_SR_expanded_six_hump_camel_back final : public Function
		{
		protected:
			virtual void initialize_(Environment *env) override;
			virtual void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_MMOP_F07 = cec2015::F7_SR_expanded_six_hump_camel_back;
}
#endif // !OFEC_F7_SR_EXPANDED_SIX_HUMP_CAMEL_BACK_H




