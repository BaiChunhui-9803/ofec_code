/********* Begin Register Information **********
{
	"name": "LIPS",
	"identifier": "LIPS",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created at 21 January 2017
// Modifed by WangJunChen at 2 Dec 2019 

/* ---------------------------------------------------------------------------------------
* @article{qu2012distance,
  volume = {17},
  number = {3},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {387--402},
  year = {2012},
  author = {Boyang Qu and Ponnuthurai Nagaratnam Suganthan and Swagatam Das},
  title = {A distance-based locally informed particle swarm model for multimodal optimization}
}
-----------------------------------------------------------------------------------------*/

#ifndef OFEC_LIPS_H
#define OFEC_LIPS_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "lips_pop.h"

namespace ofec {
	class LIPS :virtual  public Algorithm {
		OFEC_CONCRETE_INSTANCE(LIPS)
	protected:
		std::unique_ptr<SwarmLIP> m_pop;
		size_t m_pop_size;
		bool m_use_history_nearest;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_LIPS_H

