/********* Begin Register Information **********
{
	"name": "ANDE",
	"identifier": "ANDE",
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
*********************************************************************************/
// Created on 20th May, 2020 by Junchen Wang
/*
Z. Wang et al., "Automatic Niching Differential Evolution With Contour Prediction Approach for Multimodal Optimization Problems," 
in IEEE Transactions on Evolutionary Computation, vol. 24, no. 1, pp. 114-128, Feb. 2020, doi: 10.1109/TEVC.2019.2910721.
*/

#ifndef OFEC_ANDE_H
#define OFEC_ANDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "ande_pop.h"

namespace ofec {
	class ANDE :virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(ANDE)
	protected:
		std::unique_ptr<PopANDE> m_pop;
		size_t m_pop_size;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}


#endif // !OFEC_ANDE_H

