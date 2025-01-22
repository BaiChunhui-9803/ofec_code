/********* Begin Register Information **********
{ 
	"name": "MMOP_CEC2013_F09",
	"identifier": "MMOP_CEC2013_F09",
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

#ifndef OFEC_CEC13_COMPOSITION1_H
#define OFEC_CEC13_COMPOSITION1_H

#include "composition.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	namespace cec2013 {
		class Composition1 : public Composition, public MetricsMMOP {
			OFEC_CONCRETE_INSTANCE(Composition1)
		protected:
			void addInputParameters();
			void updateOptima(Environment *env) override;
			void setFunction(Environment *env) override;
			void setRotation() override;
		};
	}
	using MMOP_CEC2013_F09 = cec2013::Composition1;
}
#endif // !OFEC_CEC13_COMPOSITION1_H

