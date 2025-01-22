/********* Begin Register Information **********
{
	"name": "CrowdingDE",
	"identifier": "CrowdingDE",
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
// updated on Mar 28, 2018 by Li Zhou
// Updated on 15th Aug, 2019 by Junchen Wang

/********************************************************************************
@inproceedings{thomsen2004multimodal,
  volume = {2},
  pages = {1382--1389},
  year = {2004},
  author = {Rene Thomsen},
  organization = {IEEE},
  title = {Multimodal optimization using crowding-based differential evolution},
  booktitle = {Proceedings of the 2004 Congress on Evolutionary Computation}
}
********************************************************************************/

#ifndef OFEC_CRDE_H
#define OFEC_CRDE_H

#include "../../../../template/classic/differential_evolution/population.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class PopCrowdingDE final : public PopulationDE<> {
	public:
		PopCrowdingDE(size_t size_pop, Environment *env);
		int evolve(Environment *env, Random *rnd) override;
	protected:
		size_t nearestNeighbour(int idx, Environment *env);
	};

	class CrowdingDE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(CrowdingDE)
	protected:
		size_t m_pop_size;

		void addInputParameters();
		void run_(Environment *env) override;
	};
}
#endif // OFEC_CRDE_H
