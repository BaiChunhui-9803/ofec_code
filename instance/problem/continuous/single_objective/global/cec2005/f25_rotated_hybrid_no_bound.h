/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F25",
	"identifier": "GOP_CEC2005_F25",
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

#ifndef OFEC_F25_ROTATED_HYBRID_NO_BOUND_H
#define OFEC_F25_ROTATED_HYBRID_NO_BOUND_H

#include "../../../init_pop_bounded.h"
#include "composition.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class RotatedHybridNoBound : public InitPopBounded, public Composition, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(RotatedHybridNoBound)
		protected:
			void addInputParameters();
			void initialize_(Environment *env) override;
			void evaluateObjective(Real *x, std::vector<Real>& obj) const override;
			void setFunction(Environment *env) override;
			void setTranslation() override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F25 = cec2005::RotatedHybridNoBound;
}
#endif // !OFEC_F25_ROTATED_HYBRID_NO_BOUND_H

