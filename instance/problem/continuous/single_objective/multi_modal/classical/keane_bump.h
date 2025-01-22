/********* Begin Register Information **********
{
	"name": "Classic_Keane_Bump",
	"identifier": "KeaneBump",
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
******************************************************************************************
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*  Paper: Minimization of Keane��s Bump Function by the Repulsive Particle Swarm and the 
*			Differential Evolution Methods, Mishra, SK (2007)
*******************************************************************************************/

#ifndef OFEC_KEANE_BUMP_H
#define OFEC_KEANE_BUMP_H

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {	
	class KeaneBump : public Continuous {
		OFEC_CONCRETE_INSTANCE(KeaneBump)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_KEANE_BUMP_H