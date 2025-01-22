/********* Begin Register Information **********
{
	"name": "Classic_bent_ciga_rotated",
	"identifier": "RotatedBentCigar",
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

#ifndef OFEC_ROTATED_BENT_CIGAR_H
#define OFEC_ROTATED_BENT_CIGAR_H

#include "bent_cigar.h"

namespace ofec {
	class RotatedBentCigar : public BentCigar {
		OFEC_CONCRETE_INSTANCE(RotatedBentCigar)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
	};
}
#endif // !OFEC_ROTATED_BENT_CIGAR_H


