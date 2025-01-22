/********* Begin Register Information **********
{
	"name": "Classic_trigonometric",
	"identifier": "Trigonometric",
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
*  
*******************************************************************************************/

#ifndef OFEC_TRIGONOMETRIX_H
#define OFEC_TRIGONOMETRIX_H

#include "../../../init_pop_bounded.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class Trigonometric : public InitPopBounded, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(Trigonometric)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif //!OFEC_TRIGONOMETRIX_H