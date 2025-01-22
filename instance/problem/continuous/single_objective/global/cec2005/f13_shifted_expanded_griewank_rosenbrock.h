/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F13",
	"identifier": "GOP_CEC2005_F13",
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
*************************************************************************/


#ifndef OFEC_F13_SHIFTED_EXPANDED_GRIEWANK_ROSENBROCK_H
#define OFEC_F13_SHIFTED_EXPANDED_GRIEWANK_ROSENBROCK_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class ShiftedExpandedGriewankRosenbrock : public Function, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(ShiftedExpandedGriewankRosenbrock)
		protected:
			void addInputParameters();
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
			void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
		};
	}
	using GOP_CEC2005_F13 = cec2005::ShiftedExpandedGriewankRosenbrock;
}
#endif // ! OFEC_F13_SHIFTED_EXPANDED_GRIEWANK_ROSENBROCK_H


