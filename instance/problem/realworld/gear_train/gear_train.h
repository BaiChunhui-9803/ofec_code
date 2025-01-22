/********* Begin Register Information **********
{
	"name": "gear-train",
	"identifier": "GearTrain",
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

/**************************************************************************
design of a gear train,  which
was introduced in
E.~Sandgren, ``Nonlinear integer and discrete programming in mechanical
design,'' in \emph{the ASME Design Technology Conf.,}, 1988, pp. 95--105.
, is to optimize the gear ratio for a compound gear train that contains three
gears. It is to be designed that the gear ratio is as close as possible to
1/6.931. For each gear, the number of teeth must be between
12 and 60.
******************************************************************************/
// update Mar 29, 2018 

#ifndef OFEC_FGEAR_TRAIN_H
#define OFEC_FGEAR_TRAIN_H

#include "../../../../core/problem/continuous/continuous.h"
#include "../../single_objective/metrics_gop.h"

namespace ofec {
	class GearTrain : public Continuous, public MetricsGOP {
		OFEC_CONCRETE_INSTANCE(GearTrain)
	protected:
		void addInputParameters() {}
		virtual void initialize_(Environment* env) override;
		void updateOptima(Environment* env) override;
		void evaluateObjective(Real *x, std::vector<Real>& obj) override;

	};
}
#endif // OFEC_FGEAR_TRAIN_H