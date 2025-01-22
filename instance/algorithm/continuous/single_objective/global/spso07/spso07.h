/********* Begin Register Information **********
{
	"name": "SPSO-07",
	"identifier": "SPSO07",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* SPSO07: C. Maurice, "Standard pso 2007 (SPSO-07)" http://www.particleswarm.info/Programs.html, 2007.
*
*********************************************************************************/
// Created: 21 September 2011
// Last modified: 21 Nov 2014

#ifndef OFEC_SPSO07_H
#define OFEC_SPSO07_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/particle_swarm_optimization/swarm.h"
#include "particle07.h"

namespace ofec {
	class SPSO07 : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(SPSO07)
	protected:
		std::unique_ptr<Swarm<Particle07>> m_pop;
		Real m_weight, m_accelerator1, m_accelerator2;
		size_t m_pop_size;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		void initPop(Environment *env);
	};

}

#endif // OFEC_SPSO07_H
