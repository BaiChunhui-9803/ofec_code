/********* Begin Register Information **********
{
	"name": "Sampling-TSP-Visualization-localHnsw",
	"identifier": "SamplingTSP_Visualization_localHnsw",
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
*  Created by Diao Yiya on 2024-05-07
*
*************************************************************************/

#ifndef OFEC_SAMPLING_TSP_VISUALIZATION_LOCAL_HNSW_H
#define OFEC_SAMPLING_TSP_VISUALIZATION_LOCAL_HNSW_H

#include "../sampling_tsp_visualization.h"

#include "../../../../../../problem/combination/basin_divisioin/init_pop_basin.h"
#include "../hnsw_basin/hnsw_basin_algorithm.h"



namespace ofec {

#define CAST_TSP_HNSW_ALG(alg) dynamic_cast<SamplingTSP_Visualization_localHnsw*>(alg)

	class SamplingTSP_Visualization_localHnsw : public SamplingTSP_Visualization,
		public HnswBasinAlgorithm {
		OFEC_CONCRETE_INSTANCE(SamplingTSP_Visualization_localHnsw)
	protected:
		void initialize_(Environment* env) override;
		void run_(Environment* env) override;
	public:
		void addInputParameters();
	};
}

#endif