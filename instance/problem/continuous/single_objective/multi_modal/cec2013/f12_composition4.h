/********* Begin Register Information **********
{
	"name": "MMOP_CEC2013_F12",
	"identifier": "MMOP_CEC2013_F12",
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
*************************************************************************
*  Paper: Benchmark Functions for CEC'2013 Special Session and Competition on Niching
*  Methods for Multimodal Function Optimization.
*******************************************************************************************/

#ifndef OFEC_CEC13_COMPOSITION4_H
#define OFEC_CEC13_COMPOSITION4_H

#include "composition.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	namespace cec2013 {
		class Composition4 : public Composition, public MetricsMMOP {
			OFEC_CONCRETE_INSTANCE(Composition4)
		protected:
			void addInputParameters();
			void updateOptima(Environment *env) override;
			void setFunction(Environment *env) override;
		};
	}
	using MMOP_CEC2013_F12 = cec2013::Composition4;
}
#endif // !OFEC_CEC13_COMPOSITION4_H



