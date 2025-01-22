/********* Begin Register Information **********
{
	"name": "Classic_Rosenbrock_rotated",
	"identifier": "RotatedRosenbrock",
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

#ifndef OFEC_ROTATED_ROSENBROCK_H
#define OFEC_ROTATED_ROSENBROCK_H

#include "rosenbrock.h"

namespace ofec {
	class RotatedRosenbrock : public Rosenbrock {
		OFEC_CONCRETE_INSTANCE(RotatedRosenbrock)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
	};
}
#endif // !OFEC_ROTATED_ROSENBROCK_H

