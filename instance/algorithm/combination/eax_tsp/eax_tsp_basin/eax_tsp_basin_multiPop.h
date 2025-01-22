/********* Begin Register Information **********
{
	"name": "EAX-TSP-basin-multipop",
	"identifier": "EaxTspBasinMultiPop",
	"tags": [  "travelling salesman problem", "single-objective" ]
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
*
*  Created by Diao Yiya on 20240514
*
*-------------------------------------------------------------------------------
*

*
*************************************************************************/

#ifndef OFEC_EAX_TSP_BASIN_MULTIPOP_H
#define OFEC_EAX_TSP_BASIN_MULTIPOP_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../core/algorithm/multi_population.h"
#include "../eax_tsp_origin/environment.h"
#include "hnsw_basin.h"
#include "eax_tsp_basin_pop.h"



namespace ofec {
#define CAST_EAXTSP_BASIN_MP(alg) dynamic_cast<EaxTspBasinMultiPop*>(alg)

	class EaxTspBasinMultiPop : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(EaxTspBasinMultiPop)

		void addInputParameters() {}
	public:
		bool terminating() override;
		SolutionBase* bestSolution() {
			return m_best.get();
		}

		void udpateBestSolution();
		
	protected:

		void initialize_(Environment* env) override;
		void run_(Environment* env) override;

	protected:
		HNSWbasin m_hnswBasin;
		MultiPopulation<PopEaxTspBasin> m_pops;
		std::shared_ptr<SolutionBase> m_best;
		
		int m_maxGen = 300;
		int m_maxStagnation = 100;
	};
}

#endif // !OFEC_CMAES_H

