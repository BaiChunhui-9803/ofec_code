/********* Begin Register Information **********
{
	"name": "AMP-GL-SaDE",
	"identifier": "AMP_GL_SaDE",
	"tags": [ "contamination source identification for water distribution network" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
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
// updated Jun 10, 2019 by Li Zhou

#ifndef OFEC_AMP_CC_CSIWDN_INSTANCE_H
#define OFEC_AMP_CC_CSIWDN_INSTANCE_H

#include "amp_cc.h"
#include "../sade_population.h"
#include "../gl_population.h"
#include "../../../../../core/instance.h"

namespace ofec {
	//using AMP_GL_SaDE = AMP_CC_CSIWDN<GLPopulation, SaDEPopulation, IndCSIWDN>;

	class AMP_GL_SaDE : virtual public AMP_CC_CSIWDN<GLPopulation, SaDEPopulation, IndCSIWDN>
	{
		OFEC_CONCRETE_INSTANCE(AMP_GL_SaDE)
		void addInputParameters() {}
	};
}

#endif // !OFEC_AMP_CC_CSIWDN_INSTANCE_H