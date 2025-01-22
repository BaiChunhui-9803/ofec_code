/********* Begin Register Information **********
{
	"name": "Classic_Shekel_modified",
	"identifier": "ModifiedShekel",
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

#ifndef OFEC_MODIFIED_SHEKEL_H
#define OFEC_MODIFIED_SHEKEL_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class ModifiedShekel : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(ModifiedShekel)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;

		Real m_a[8][5];
		Real m_c[8];
	};
}
#endif // !OFEC_MODIFIED_SHEKEL_H