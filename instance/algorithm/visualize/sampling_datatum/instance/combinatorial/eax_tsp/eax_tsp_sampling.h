/********* Begin Register Information **********
{
	"name": "EAX-TSP-data-sampling",
	"identifier": "EAX_TSP_Data_Sampling",
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
*  Created by Diao Yiya on 2024-09-16
*
*************************************************************************/

#ifndef OFEC_SAMPLING_EAX_TSP_INSTANCE_H
#define OFEC_SAMPLING_EAX_TSP_INSTANCE_H


#include "../../../sampling_algorithm_multipop.h"
#include "../../../../../../algorithm/combination/eax_tsp/eax_tsp_origin/eax_tsp_alg.h"
namespace ofec {
	//using EAX_TSP_Data_Sampling = SamplingAlgorithm;
	class EAX_TSP_Data_Sampling : public SamplingAlgorithm<EAX_TSP> {
		OFEC_CONCRETE_INSTANCE(EAX_TSP_Data_Sampling)
	};
}


#endif