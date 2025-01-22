/********* Begin Register Information **********
{
	"name": "TSC2-sampling",
	"identifier": "TSC2_Sampling",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Created by Diao Yiya on 2024-10-03
*
*************************************************************************/

#ifndef OFEC_SAMPLING_TSC2_INSTANCE_H
#define OFEC_SAMPLING_TSC2_INSTANCE_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/multi_modal/tsc2/tsc2.h"
namespace ofec {

	class TSC2_Sampling : public SamplingAlgorithm<TSC2> {
		OFEC_CONCRETE_INSTANCE(TSC2_Sampling)
	};
}
#endif

