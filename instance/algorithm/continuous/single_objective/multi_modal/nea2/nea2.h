/********* Begin Register Information **********
{
	"name": "NEA2",
	"identifier": "NEA2",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Implemented by Junchen Wang (wangjunchen.chris@gmail.com) at 2020/9/23
*********************************************************************************/

/*************************************************
@inproceedings{preuss2012improved,
  pages = {386--395},
  year = {2012},
  author = {Mike Preuss},
  organization = {Springer},
  title = {Improved topological niching for real-valued global optimization},
  booktitle = {European Conference on the Applications of Evolutionary Computation}
}
*************************************************/

#ifndef OFEC_NEA2_H
#define OFEC_NEA2_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../global/cma_es/cma_es_pop.h"

namespace ofec {
	class NEA2 : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(NEA2)
	protected:
		MultiPopulation<PopCMA_ES> m_subpops;
		size_t m_pop_size;
		bool m_without_restart;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		void addSubpops(size_t num_samples, Environment *env);
		bool stopTolFun(PopCMA_ES &subpop);
	};
}

#endif // !OFEC_NEA2_H
