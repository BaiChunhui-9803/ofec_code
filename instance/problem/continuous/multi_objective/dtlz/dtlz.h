/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 December 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)

/***********************************************************************************
Deb, K., Thiele, L., Laumanns, M., & Zitzler, E. (2005). 
Scalable test problems for evolutionary multiobjective optimization. 
In Evolutionary multiobjective optimization (pp. 105-145). Springer London.
************************************************************************************/

#ifndef OFEC_DTLZ_H
#define OFEC_DTLZ_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	typedef std::vector<Real> TObjVec;
	typedef std::vector<TObjVec> TFront;

	class DTLZ : public Continuous, public MetricsMOP  {
	protected:
		void initialize_();
		void generateParetoFront();
		void loadParetoFront();
		void generateRecursive(TFront *pf, TObjVec *pt, size_t number_objectives, size_t left, size_t total, size_t element);
		void generateWeight(TFront *pf, size_t M, size_t p);
		void generateOnelayerPF(std::ostream &os, const std::string &problem_name, int M, int p);
		void generateTwolayersPF(std::ostream &os, const std::string &problem_name, int M, int outside_p, int inside_p);
	};
}

#endif //OFEC_DTLZ_H