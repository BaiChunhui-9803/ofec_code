/********* Begin Register Information **********
{
	"name": "Classic_six_hump_camel_back",
	"identifier": "SixHumpCamelBack",
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
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
//Stoean, C.; Preuss, M.; Stoean, R.; Dumitrescu, D., 
// "Multimodal Optimization by Means of a Topological Species Conservation Algorithm," 
// Evolutionary Computation, IEEE Transactions on , vol.14, no.6, pp.842,864, Dec. 2010
//doi: 10.1109/TEVC.2010.204166
#ifndef OFEC_SIX_HUMP_H
#define OFEC_SIX_HUMP_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class SixHumpCamelBack : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(SixHumpCamelBack)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_SIX_HUMP_H