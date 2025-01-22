/********* Begin Register Information **********
{
	"name": "Classic_discus_rotated",
	"identifier": "RotatedDiscus",
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

#ifndef OFEC_ROTATED_DISCUS_H
#define OFEC_ROTATED_DISCUS_H

#include "discus.h"

namespace ofec {
	class RotatedDiscus : public Discus {
		OFEC_CONCRETE_INSTANCE(RotatedDiscus)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
	};
}
#endif // !OFEC_ROTATED_DISCUS_H


