/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Junchen Wang
* Email: changhe.lw@gmail.com, wangjunchen.chris@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************
*  Paper: Benchmark Functions for CEC��2013 Special Session and Competition on Niching
*  Methods for Multimodal Function Optimization.
*******************************************************************************************/

#ifndef OFEC_CEC13_HIMMENBLAU_H
#define OFEC_CEC13_HIMMENBLAU_H

#include "../classical/himmenblau.h"

namespace ofec {
	namespace cec2013 {
		// An inverted version of Himmelblau function
		class Himmenblau : public ofec::Himmenblau {
			OFEC_CONCRETE_INSTANCE(Himmenblau)
		protected:
			void addInputParameters() {}
			void initialize_(Environment *env) override;
			void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
		};
	}
}
#endif // !OFEC_CEC13_F4_HIMMENBLAU_H

