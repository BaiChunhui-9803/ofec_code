/********* Begin Register Information **********
{
	"name": "EAX-TSP-basin",
	"identifier": "EaxTspBasin",
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
*  Created by Diao Yiya on 2023 11 13
*
*-------------------------------------------------------------------------------
*

*
*************************************************************************/

#ifndef OFEC_EAX_TSP_BASIN_H
#define OFEC_EAX_TSP_BASIN_H

#include "../../../../../core/algorithm/algorithm.h"
#include "hnsw_basin.h"
#include "eax_tsp_basin_pop.h"


namespace ofec {
	class EaxTspBasin : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(EaxTspBasin)

		void addInputParameters();
	public:
		bool terminating() override {
			return Algorithm::terminating() || m_pop->terminationCondition();
			//return Algorithm::terminating() || m_pop->terminationCondition();
		}
	protected:
		std::unique_ptr<PopEaxTspBasin> m_pop;

		void initialize_(Environment* env) override;
		void run_(Environment* env) override;

	protected:
		HNSWbasin m_hnswBasin;
		
	};
}

#endif // !OFEC_CMAES_H

