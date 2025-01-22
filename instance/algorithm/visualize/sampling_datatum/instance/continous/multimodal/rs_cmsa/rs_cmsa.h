/********* Begin Register Information **********
{
	"name": "RS_CMSA-sampling",
	"identifier": "RS_CMSA_Sampling",
	"tags": [ "continuous", "single-objective" ],
	"dependency on libraries": [ "Eigen" ]
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

#ifndef OFEC_SAMPLING_RS_CMSA_INSTANCE_H
#define OFEC_SAMPLING_RS_CMSA_INSTANCE_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/multi_modal/rs_cmsa/rs_cmsa.h"
namespace ofec {

	class RS_CMSA_Sampling : public SamplingAlgorithm<RS_CMSA> {
		OFEC_CONCRETE_INSTANCE(RS_CMSA_Sampling)
	};

}
#endif

