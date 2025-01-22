/********* Begin Register Information **********
{
	"name": "Classic_Szu",
	"identifier": "Szu",
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

#ifndef OFEC_SZU_H
#define OFEC_SZU_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class Szu : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(Szu)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_SZU_H