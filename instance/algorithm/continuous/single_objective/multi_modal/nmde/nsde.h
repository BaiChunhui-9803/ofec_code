/********* Begin Register Information **********
{
	"name": "NSDE",
	"identifier": "NSDE",
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
// updated Mar 28, 2018 by Li Zhou

/******************************************************************************************
@article{qu2012differential,
  volume = {16},
  number = {5},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {601--614},
  year = {2012},
  author = {Boyang Qu and Ponnuthurai Nagaratnam Suganthan and Jing Liang},
  title = {Differential evolution with neighborhood mutation for multimodal optimization}
}
******************************************************************************************/

#ifndef OFEC_NSDE_H
#define OFEC_NSDE_H

#include "../../../../template/classic/differential_evolution/population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include <list>

namespace ofec {
	class PopNSDE final : public PopulationDE<>	{
	public:
		PopNSDE(size_t size_pop, size_t cluster_size, Environment *env);
		void selectSubpop(Environment *env);
		int evolve(Environment *env, Random *rnd) override;
		const std::vector<std::list<std::pair<Real, int>>>& dist() const { return m_dis; }
	protected:
		int m_m;                                //size of neighborhood
		std::vector<std::list<std::pair<Real, int>>> m_dis;  //save Solutions' distance
		std::vector<int> m_seed;                     //best fittness of neighborhood
		std::vector<int> m_order_list;
	};

	class NSDE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(NSDE)
	protected:
		std::unique_ptr<PopNSDE> m_pop;
		size_t m_pop_size, m_cluster_size;

		void addInputParameters();
		void run_(Environment *env) override;
		void initPop(size_t pop_size, Environment *env);
	};
}


#endif // NSDE_H
