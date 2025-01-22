/********* Begin Register Information **********
{
	"name": "MOP_DTLZ3",
	"identifier": "DTLZ3",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 12 JAN 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)

#ifndef OFEC_DTLZ3_H
#define OFEC_DTLZ3_H

#include "dtlz.h"

namespace ofec {
	class DTLZ3 : public DTLZ {
	protected:
		void evaluateObjective(Real *x, std::vector<Real> &obj) override;
	};
}

#endif //OFEC_DTLZ3_H  