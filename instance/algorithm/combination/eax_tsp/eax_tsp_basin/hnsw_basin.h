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

#ifndef OFEC_HNSW_BASIN_H
#define OFEC_HNSW_BASIN_H
#include "../../../../../utility/hnsw/hnsw_nbn/hnsw_model.h"

#include "../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class HNSWbasin {
		nbn::HnswModel m_model;
		std::vector<int> m_belongBasinId;
		std::vector<std::vector<int>> m_basinIds;
		int m_numNeighbors = 20;
	public:
		const nbn::HnswModel& hnswModel()const {
			return m_model;
		}
		const std::vector<int>& belongBasinId()const{
			return m_belongBasinId;
		}
		const std::vector<std::vector<int>>& basinIds()const {
			return m_basinIds;
		}


		void initialize(const nbn::HnswModel& model,
			const std::vector<int>& belongBasinId,
			const std::vector<std::vector<int>>& basinIds);
		// threadSafe
		int calculateBasinId(SolutionBase& sol)const;
		
	};
}


#endif