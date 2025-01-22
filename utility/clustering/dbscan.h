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
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
* Ester M, Kriegel H P, Sander J, et al. 
* Density-based spatial clustering of applications with noise
* Int. Conf. knowledge discovery and data mining. 1996, 240(6).
* source code is from: https://github.com/james-yoo/DBSCAN
* The algorithm had implemented with pseudocode described in [wiki](https://en.wikipedia.org/wiki/DBSCAN), 
* but it is not optimised. 
*-------------------------------------------------------------------------------
* Constructed by Qingshan Tan at 2023/4/6
*********************************************************************************/

#ifndef OFEC_DBSCAN_H
#define OEFC_DBSCAN_H

#include <algorithm>
#include <vector>
#include <set>
#include <numeric>

#include"../../core/algorithm/population.h"
#include"../../core/problem/solution.h"

#include<fstream>

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

namespace ofec {

	typedef struct Point_
	{
		std::vector<Real> m_variables;  // point position
		int clusterID;  // clustered ID
	}Point;

	class DBSCAN {
	public:
		std::vector<std::unique_ptr<Point>> m_points;
	private:
		unsigned int m_pointSize;
		unsigned int m_minPoints;
		float m_epsilon;
	public:
		DBSCAN(unsigned int minPts, float eps, const std::vector<std::vector<Real>> &points) {
			m_minPoints = minPts;
			m_epsilon = eps;
			m_pointSize = points.size();
			for (size_t i = 0; i < points.size(); ++i) {
				m_points.emplace_back(new Point);
			}
			for (size_t i = 0; i < points.size(); ++i) {
				m_points[i]->m_variables=points[i];
				m_points[i]->clusterID = UNCLASSIFIED;
			}
		}
		~DBSCAN() {}

		int run();
		std::vector<int> calculateCluster(const Point &point);
		int expandCluster(Point &point, int clusterID);
		inline Real calculateDistance(const Point &pointCore, const Point &pointTarget);

		int getTotalPointSize() { return m_pointSize; }
		int getMinimumClusterSize() { return m_minPoints; }
		int getEpsilonSize() { return m_epsilon; }
	};
}

#endif // !OFEC_DBSCAN_H