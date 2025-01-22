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

#ifndef OFEC_HNSW_BASIN_ALGORITHM_H
#define OFEC_HNSW_BASIN_ALGORITHM_H

#include "../../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../../../../../../../utility/hnsw/hnsw_nbn/hnsw_nbn.h"

namespace ofec {

#define CAST_HNSW_B(alg) dynamic_cast<HnswBasinAlgorithm*>(alg)

	class HnswBasinAlgorithm : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(HnswBasinAlgorithm)

	protected:


		bool m_localNBN = false;
		nbn::HnswModel m_model;
		std::vector<Real> m_fitness;
		std::vector<int> m_belong;
		std::vector<Real> m_dis2parent;



	protected:
	//	void initialize_(Environment* env) override {}
		void run_(Environment* env) override {}

	public:
		void addInputParameters() {}


		const nbn::HnswModel& hnswModel() const {
			return m_model;
		}
		nbn::HnswModel& hnswModel() {
			return m_model;
		}
		const std::vector<ofec::Real>& dis2parent()const {
			return m_dis2parent;
		}
		const std::vector<int>& belong()const {
			return m_belong;
		}


	};
}

#endif