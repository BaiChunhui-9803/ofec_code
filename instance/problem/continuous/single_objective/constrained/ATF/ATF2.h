/********* Begin Register Information **********
{
	"name": "ATF2",
	"identifier": "ATF2",
	"tags": [ "single-objective", "continuous", "constrained" ]
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


#ifndef OFEC_ATF2_H
#define OFEC_ATF2_H

#include "../cop_base.h"

namespace ofec {
	namespace ATF {
		class _2 final: public CopBase
		{
			OFEC_CONCRETE_INSTANCE(_2)
		public:
			void initialize_(Environment* env) override;
		protected:
			void evaluateObjectiveAndConstraint(Real* x, std::vector<Real>& obj, std::vector<Real>& con) override;
			void addInputParameters();

		};
	}
	using ATF2 = ATF::_2;
}
#endif // ! ATF2
