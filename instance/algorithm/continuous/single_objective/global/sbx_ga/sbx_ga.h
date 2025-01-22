/********* Begin Register Information **********
{
	"name": "SBX-GA",
	"identifier": "SBX_GA",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com
* Language: C++
* -----------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

/*************************************************************************
* Deb K, Agrawal R B. 
* Simulated binary crossover for continuous search space[J]. 
* Complex systems, 1995, 9(2): 115-148.
* -----------------------------------------------------------------------
* Created: 29 April 2022
*************************************************************************/

#ifndef OFEC_SBX_GA_H
#define OFEC_SBX_GA_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/genetic_algorithm/sbx_pop.h"

namespace ofec {
	class SBX_GA : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(SBX_GA)
	protected:
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		PopSBX<> m_pop;
		Population<Solution<>> m_pop_and_off;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_SBX_GA_H
