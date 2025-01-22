/********* Begin Register Information **********
{
	"name": "EMO_MMO_sampling",
	"identifier": "EMO_MMO_Sampling",
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

#ifndef OFEC_SAMPLING_EMO_MMO_INSTANCE_H
#define OFEC_SAMPLING_EMO_MMO_INSTANCE_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/multi_modal/emo_mmo/emo_mmo.h"
namespace ofec {

	class EMO_MMO_Sampling : public SamplingAlgorithm<EMO_MMO> {
		OFEC_CONCRETE_INSTANCE(EMO_MMO_Sampling)
	};
}
#endif
