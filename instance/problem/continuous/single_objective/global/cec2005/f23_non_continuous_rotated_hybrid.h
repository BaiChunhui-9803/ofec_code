/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F23",
	"identifier": "GOP_CEC2005_F23",
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

#ifndef OFEC_F23_NON_CONTINUOUS_ROTATED_HYBRID_H
#define OFEC_F23_NON_CONTINUOUS_ROTATED_HYBRID_H

#include "composition.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class NonContinuousRotatedHybrid : public Composition, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(NonContinuousRotatedHybrid)
		protected:
			void addInputParameters();
			void evaluateObjective(Real *x, std::vector<Real>& obj) const override;
			void setFunction(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F23 = cec2005::NonContinuousRotatedHybrid;
}
#endif // !OFEC_F23_NON_CONTINUOUS_ROTATED_HYBRID_H
