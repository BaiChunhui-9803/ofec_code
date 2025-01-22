/********* Begin Register Information **********
{
	"name": "Classic_elliptic",
	"identifier": "Elliptic",
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

#ifndef OFEC_ELLIPTIC_H
#define OFEC_ELLIPTIC_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {	
	class Elliptic : public Function, public MetricsGOP {
		OFEC_CONCRETE_INSTANCE(Elliptic)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
	};	
}

#endif // !OFEC_ELLIPTIC_H
