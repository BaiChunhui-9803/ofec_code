/********* Begin Register Information **********
{
	"name": "MMOP_CEC2015_F04",
	"identifier": "CEC2015_MMOP_F04",
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


#ifndef OFEC_F4_SR_EXPANDED_DECREASING_MINIMA_H
#define OFEC_F4_SR_EXPANDED_DECREASING_MINIMA_H

#include "function.h"

namespace ofec {
	namespace cec2015 {
		class F4_SR_expanded_decreasing_minima final : public Function
		{
		public:
			//F4_SR_expanded_decreasing_minima(const ParameterMap &v);
			//F4_SR_expanded_decreasing_minima(const std::string &name, size_t size_var, size_t size_obj);
		protected:
			virtual void initialize_(Environment *env) override;
			virtual void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_MMOP_F04 = cec2015::F4_SR_expanded_decreasing_minima;
}
#endif // !OFEC_F4_SR_EXPANDED_DECREASING_MINIMA_H

