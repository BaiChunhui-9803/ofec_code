/********* Begin Register Information **********
{
	"dependency on libraries": [ "Eigen" ],
	"name": "EA4eig-sampling",
	"identifier": "EA4eig_Sampling",
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

#ifndef OFEC_SAMPLING_EA4EIG_H
#define OFEC_SAMPLING_EA4EIG_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/global/ea4eig/ea4eig.h"

namespace ofec {
	class EA4eig_Sampling : public SamplingAlgorithm<EA4eig> {
		OFEC_CONCRETE_INSTANCE(EA4eig_Sampling)
	};
}
#endif

