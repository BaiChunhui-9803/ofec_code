/********* Begin Register Information **********
{
	"name": "GOP_CEC2015_F08",
	"identifier": "CEC2015_GOP_F08",
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
*************************************************************************
*  Paper : Problem Definitions and Evaluation Criteria for the CEC2015
*  Competition on Learning-based Real-Parameter Single Objective
*  Optimization.
************************************************************************/
#ifndef OFEC_F8_HYBRID3_H
#define OFEC_F8_HYBRID3_H

#include "hybrid.h"

namespace ofec {
	namespace cec2015 {
		class F8_hybrid3 final : public Hybrid {
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction() override;
		};
	}
	using CEC2015_GOP_F08 = cec2015::F8_hybrid3;
}

#endif // !OFEC_F8_HYBRID3_H








