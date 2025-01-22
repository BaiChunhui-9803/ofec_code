/********* Begin Register Information **********
{
	"name": "DBLMA",
	"identifier": "DBLMA",
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
*********************************************************************************/

/******************************************************************************
@article{segura2015novel,
  volume = {46},
  number = {12},
  journal = {IEEE Transactions on Cybernetics},
  pages = {3233--3246},
  year = {2015},
  author = {Carlos Segura and Carlos A. Coello Coello and Eduardo Segredo and Arturo Hernandez Aguirre},
  title = {A novel diversity-based replacement strategy for evolutionary algorithms}
}
*-------------------------------------------------------------------------------
* Created on 6th April, 2022 by Junchen Wang
*********************************************************************************/

#ifndef OFEC_DBLMA_H
#define OFEC_DBLMA_H

#include "../../../../template/classic/genetic_algorithm/sbx_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	/* Diversity-Based Lamarckian Memetic Algorithm */
	class DBLMA : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(DBLMA)
	protected:
		enum class SurvivorSelectionScheme { kMulti, kMultiDynamic };

		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		SurvivorSelectionScheme m_sss;
		PopSBX<> m_pop;
		Population<Solution<>> m_pop_and_off;
		std::vector<Real> m_step_size;
		size_t m_max_iter_stagnant;
		std::unique_ptr<Solution<>> m_neighbor;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		void localSearch(Population<Solution<>>& pop, Environment *env);
	};
}

#endif // !OFEC_DBLMA_H
