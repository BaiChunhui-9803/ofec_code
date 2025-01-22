/********* Begin Register Information **********
[
	{ "name": "MOEA_F1", "identifier": "MOEA_F1", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F2", "identifier": "MOEA_F2", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F3", "identifier": "MOEA_F3", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F4", "identifier": "MOEA_F4", "problem tags": [ "MOP", "ConOP" ] }, 
	{ "name": "MOEA_F5", "identifier": "MOEA_F5", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F6", "identifier": "MOEA_F6", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F7", "identifier": "MOEA_F7", "problem tags": [ "MOP", "ConOP" ] }, 
	{ "name": "MOEA_F8", "identifier": "MOEA_F8", "problem tags": [ "MOP", "ConOP" ] },
	{ "name": "MOEA_F9", "identifier": "MOEA_F9", "problem tags": [ "MOP", "ConOP" ] }
]
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 December 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)

/*********************************************************************************************************************
  H. Li and Q. Zhang
  Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization Problems with Complicated Pareto Sets
  Technical Report CES-476, Department of Computer Science, University of Essex, 2007
**********************************************************************************************************************/

#ifndef OFEC_MOEA_F_H
#define OFEC_MOEA_F_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class MOEA_F : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) override;

		int getDtype() const { return m_dtype; }
		int getPtype() const { return m_ptype; }
		int getLtype() const { return m_ltype; }
		// control the PF shape
		void alphaFunction(Real alpha[], const Real *x, int dim, int type);
		// control the distance
		Real betaFunction(const std::vector<Real> & x, int type);
		// control the PS shape of 2-d instances
		Real paretoSetFunc2(const Real &x, const Real &t1, size_t dim, int type, int css);
		// control the PS shapes of 3-D instances
		Real paretoSetFunc3(const Real &x, const Real &t1, const Real &t2, int dim, int type);
		void LoadParetoFront();

		int m_dtype, m_ptype, m_ltype;
	};

	using MOEA_F1 = MOEA_F;
	using MOEA_F2 = MOEA_F;
	using MOEA_F3 = MOEA_F;
	using MOEA_F4 = MOEA_F;
	using MOEA_F5 = MOEA_F;
	using MOEA_F6 = MOEA_F;
	using MOEA_F7 = MOEA_F;
	using MOEA_F8 = MOEA_F;
	using MOEA_F9 = MOEA_F;
}
#endif //MOEA_F_H