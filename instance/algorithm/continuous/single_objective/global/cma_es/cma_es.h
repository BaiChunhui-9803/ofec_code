/********* Begin Register Information **********
{
	"name": "CMA-ES",
	"identifier": "CMA_ES",
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
*  Created by Junchen Wang on Oct. 29, 2018.
*************************************************************************/

/********************************************************************
@article{hansen2001completely,
  author = {Nikolaus Hansen and Andreas Ostermeier},
  journal = {Evolutionary Computation},
  title = {Completely derandomized self-adaptation in evolution strategies},
  year = {2001},
  month = {06},
  number = {2},
  pages = {159-195},
  volume = {9}
}
********************************************************************/

#ifndef OFEC_CMAES_H
#define OFEC_CMAES_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class CMA_ES : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(CMA_ES)
	protected:
		size_t m_pop_size;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_CMAES_H

