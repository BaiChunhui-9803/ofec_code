/********* Begin Register Information **********
{
	"name": "NL_SHADE_RSP_MID-sampling",
	"identifier": "NL_SHADE_RSP_MID_Sampling",
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

#ifndef OFEC_SAMPLING_NL_SHADE_RSP_MID_H
#define OFEC_SAMPLING_NL_SHADE_RSP_MID_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/global/nl_shade_rsp_mid/nl_shade_rsp_mid.h"

namespace ofec {
	class NL_SHADE_RSP_MID_Sampling : public SamplingAlgorithm<NL_SHADE_RSP_MID> {
		OFEC_CONCRETE_INSTANCE(NL_SHADE_RSP_MID_Sampling)
	};
}
#endif

