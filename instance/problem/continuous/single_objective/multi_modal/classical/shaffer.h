/********* Begin Register Information **********
{
	"name": "Classic_Shaffer",
	"identifier": "Shaffer",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/

#ifndef OFEC_SHAFFER_H
#define OFEC_SHAFFER_H

#include "../../function.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class Shaffer : public Function, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(Shaffer)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real> &obj) override;
	};
}
#endif // !OFEC_SHAFFER_H