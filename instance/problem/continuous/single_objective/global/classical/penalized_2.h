/********* Begin Register Information **********
{
	"name": "Classic_penalized_2",
	"identifier": "Penalized_2",
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
*************************************************************************/

#ifndef OFEC_PENALIZED_2_H
#define OFEC_PENALIZED_2_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	class Penalized_2 : public Function, public MetricsGOP {
		OFEC_CONCRETE_INSTANCE(Penalized_2)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real> &obj) const override;
		Real u(Real x, Real a, Real k, Real m)const;
	};
}
#endif // !OFEC_PENALIZED_2_H
