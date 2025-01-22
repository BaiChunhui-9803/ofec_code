/********* Begin Register Information **********
{
	"name": "MOP_DTLZ7",
	"identifier": "DTLZ7",
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

//Note: It was named after DTLZ7 in [Deb2001] and source file "dtlz_c_source.tar.gz"

// Created: 7 Oct 2018 by Junchen Wang (wangjunchen.chris@gmail.com)
// Modified: 

#ifndef OFEC_DTLZ7_H
#define OFEC_DTLZ7_H

#include "dtlz.h"

namespace ofec {
	class DTLZ7 : public DTLZ {
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif // !OFEC_DTLZ7_H
