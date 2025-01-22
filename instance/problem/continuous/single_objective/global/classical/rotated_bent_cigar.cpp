
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

#include "rotated_bent_cigar.h"

namespace ofec {
	void RotatedBentCigar::initialize_(Environment *env) {
		BentCigar::initialize_(env);
		loadRotation("/instance/problem/continuous/single_objective/global/classical/");
	}
}