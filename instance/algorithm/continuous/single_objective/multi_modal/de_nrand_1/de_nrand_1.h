/********* Begin Register Information **********
{
	"name": "DE/nrand/1",
	"identifier": "DE_nrand_1",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
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

/********************************************************
* @inproceedings{epitropakis2011finding,
  pages = {1--8},
  year = {2011},
  author = {Michael G. Epitropakis and Vassilis P. Plagianakos and Michael N. Vrahatis},
  organization = {IEEE},
  title = {Finding multiple global optima exploiting differential evolution's niching capability},
  booktitle = {2011 IEEE Symposium on Differential Evolution (SDE)}
}
**********************************************************/

#ifndef OFEC_DE_NRAND_1_H
#define OFEC_DE_NRAND_1_H

#include "de_nrand_1_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class DE_nrand_1 :virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(DE_nrand_1)
	protected:
		std::unique_ptr<PopDE_nrand_1> m_pop;
		size_t m_pop_size;
		Real m_scaling_factor, m_crossover_rate;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		virtual void initPop(Environment *env);
	};
}

#endif // !OFEC_DE_NRAND_1_H

