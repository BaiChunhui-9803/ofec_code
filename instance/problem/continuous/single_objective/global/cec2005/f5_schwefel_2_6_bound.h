/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F05",
	"identifier": "GOP_CEC2005_F05",
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


#ifndef OFEC_F5_SCHWEFEL_2_6_BOUND_H
#define OFEC_F5_SCHWEFEL_2_6_BOUND_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class Schwefel_2_6_Bound : public Function, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(Schwefel_2_6_Bound)
		protected:
			void addInputParameters();
 			void initialize_(Environment *env) override;
			void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
			void loadData(const std::string &path);

			std::vector<Real> m_b;
			std::vector<std::vector<Real>> m_a;
		};
	}
	using GOP_CEC2005_F05 = cec2005::Schwefel_2_6_Bound;
}
#endif // ! OFEC_F5_SCHWEFEL_2_6_BOUND_H
