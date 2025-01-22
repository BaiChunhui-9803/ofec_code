
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

#include "rotated_rosenbrock.h"

namespace ofec {
	void RotatedRosenbrock::initialize_(Environment *env) {
		Rosenbrock::initialize_(env);
		setConditionNumber(1);
		loadRotation("/instance/problem/continuous/single_objective/global/classical/");	
	}
}