/********* Begin Register Information **********
{
	"name": "Sampling-TSP-Visualization",
	"identifier": "SamplingTSP_Visualization",
	"tags": [ "travelling salesman problem" ]
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
*  Created by Diao Yiya on 2024-04-02
*
*************************************************************************/

#ifndef OFEC_SAMPLING_TSP_VISUALIZATION_H
#define OFEC_SAMPLING_TSP_VISUALIZATION_H


#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class SamplingTSP_Visualization : virtual public Algorithm{
		OFEC_CONCRETE_INSTANCE(SamplingTSP_Visualization)

	protected:


		std::vector<std::vector<std::vector<int>>> m_solIds;
		std::string m_file_dir;
		int m_numRun = 0;
		int m_numIter = 0;
	protected:
		void initialize_(Environment* env) override;
		void run_(Environment* env) override;

	public:
		void addInputParameters();
		
		
	};
}

#endif