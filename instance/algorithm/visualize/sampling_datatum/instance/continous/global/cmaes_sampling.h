/********* Begin Register Information **********
{
	"name": "CMEAS-data-sampling",
	"identifier": "CMEAS_Data_Sampling",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Created by Diao Yiya on 2024-09-016
*
*************************************************************************/

#ifndef OFEC_SAMPLING_CMAES_INSTANCE_H
#define OFEC_SAMPLING_CMAES_INSTANCE_H


#include "../../../sampling_algorithm_multipop.h"
#include "../../../../../../algorithm/continuous/single_objective/global/cma_es/cma_es.h"
namespace ofec {
	using CMEAS_Data_Sampling = SamplingAlgorithm<CMA_ES>;

	//using CMEAS_Data_Sampling = SamplingAlgorithm;
}


#endif