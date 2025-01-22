/********* Begin Register Information **********
{
	"name": "ADOPT-GL-SaDE",
	"identifier": "ADOPT_GL_SaDE",
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

#ifndef OFEC_ADOPT_COMBINATION_H
#define OFEC_ADOPT_COMBINATION_H

#include "adopt.h"
#include "../sade_population.h"
#include "../gl_population.h"

namespace ofec {

	using ADOPT_GL_SaDE = ADOPT<GLPopulation, SaDEPopulation, IndCSIWDN>;

	//class ADOPT_combination : public ADOPT<GLPopulation, SaDEPopulation, IndCSIWDN>
	//{
	//public:
	//	ADOPT_combination(param_map & v) : ADOPT(v) {}
	//};

}

#endif