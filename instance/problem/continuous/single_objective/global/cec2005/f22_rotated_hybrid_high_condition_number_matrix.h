/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F22",
	"identifier": "GOP_CEC2005_F22",
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

#ifndef OFEC_F22_ROTATED_HYBRID_HIGH_CONDITION_NUMBER_MATRIX_H
#define OFEC_F22_ROTATED_HYBRID_HIGH_CONDITION_NUMBER_MATRIX_H

#include "composition.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class RotatedHybridHighConditionNumberMatrix : public Composition, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(RotatedHybridHighConditionNumberMatrix)
		protected:
			void addInputParameters();
			void evaluateObjective(Real *x, std::vector<Real>& obj) const override;
			void setFunction(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F22 = cec2005::RotatedHybridHighConditionNumberMatrix;
}
#endif // !OFEC_F22_ROTATED_HYBRID_HIGH_CONDITION_NUMBER_MATRIX_H
