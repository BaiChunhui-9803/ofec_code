/********* Begin Register Information **********
{
	"name": "Classic_Schwefel_modified",
	"identifier": "ModifiedSchwefel",
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

#ifndef OFEC_MODIFIED_SCHWEFEL_H
#define OFEC_MODIFIED_SCHWEFEL_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	class ModifiedSchwefel : public Function, public MetricsGOP {
		OFEC_CONCRETE_INSTANCE(ModifiedSchwefel)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
	};

}
#endif // !OFEC_MODIFIED_SCHWEFEL_H
