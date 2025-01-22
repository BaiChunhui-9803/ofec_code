/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Xiaofang Wu
* Email: changhe.lw@google.com Or wuxiaofang@cug.edu.cn
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
/*********************************************************************************************************************
  S. Huband, L. Barone, L. While, and P. Hingston
  A scalable multi-objective test problem toolkit. Lecture Notes inComputer Science,2005, 3410:280-295.
  **********************************************************************************************************************/
// Created: 5 August 2019


#ifndef WFG_H
#define WFG_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class WFG : public Continuous, public MetricsMOP  {
	protected:
		int m_k;				//position parameter(k) should be a multiple of M - 1, distance parameter(l) is N-k
		Real m_D;				//distance parameter
		std::vector<Real> m_S;	//shape parameter(scaling constants)(M)
		std::vector<Real> m_A;	//degeneracy constants(M-1)
		std::vector<Real> m_h;	//shape function
		
	protected:	
		void initialize_() override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) override;

		void set_k(int k) { m_k = k; }
		void loadParetoFront();

		std::vector<Real> calculateX(const std::vector<Real> &t_p);
		std::vector<Real> normalise(Real *z);

		virtual void t1(std::vector<Real> &y) {}
		virtual void t2(std::vector<Real> &y) {}
		virtual void t3(std::vector<Real> &y) {}
		virtual void t4(std::vector<Real> &y) {}
		virtual void shape(std::vector<Real> &y) {}

		Real linear(const std::vector<Real> &x, int m);	//The linear shape function. (m is indexed from 1.)
		Real convex(const std::vector<Real> &x, int m);	//The convex shape function. (m is indexed from 1.)		
		Real concave(const std::vector<Real> &x, int m);	//The concave shape function. (m is indexed from 1.)		
		Real mixed(const std::vector<Real> &x, int A, Real alpha);	//The mixed convex/concave shape function		
		Real disc(const std::vector<Real> &x, int A, Real alpha, Real beta);	//The disconnected shape function
		Real correctTo01(const Real &a, Real epsilon = 1.0e-10);

		std::vector<Real> subvector(const std::vector<Real> &v, int head, int tail);

		Real bPoly(Real y, Real alpha);
		Real bFlat(Real y, Real A, Real B, Real C);
		Real bParam(Real y, Real u, Real A, Real B, Real C);
		Real sLinear(Real y, Real A);
		Real sDecept(Real y, Real A, Real B, Real C);
		Real sMulti(Real y, int A, Real B, Real C);
		Real rSum(const std::vector<Real> &y, const std::vector<Real> &w);
		Real rNonsep(const std::vector<Real> &y, int A);

	};
}
#endif //MOEA_FBASE_H