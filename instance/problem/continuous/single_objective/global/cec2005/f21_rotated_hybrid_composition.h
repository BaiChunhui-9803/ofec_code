/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F21",
	"identifier": "GOP_CEC2005_F21",
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

#ifndef OFEC_F21_ROTATED_HYBRID_COMPOSITION_H
#define OFEC_F21_ROTATED_HYBRID_COMPOSITION_H

#include "composition.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class RotatedHybridComposition3 : public Composition, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(RotatedHybridComposition3)
		protected:
			void addInputParameters();
			void evaluateObjective(Real *x, std::vector<Real>& obj) const override;
			void setFunction(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F21 = cec2005::RotatedHybridComposition3;
}
#endif // !OFEC_F21_ROTATED_HYBRID_COMPOSITION_H
