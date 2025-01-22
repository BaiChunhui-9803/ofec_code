/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#include "dynde_Solution.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include <cmath>

namespace ofec {
	int IndDynDE::brownian(const SolutionType &best, Real sigma, Problem *pro, Algorithm *alg, Random *rnd) {
		for (size_t i = 0; i < variable().size(); ++i) {
			variable()[i] = (best.variable()[i]) + rnd->normal.nextNonStd(0, sigma);
		}
		CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
		return evaluate(pro, alg);
	}

	int IndDynDE::quantum(const SolutionType &best, Real rcloud, Problem *pro, Algorithm *alg, Random *rnd) {
		const auto dim = best.variable().size();
		std::vector<Real> x(dim, 0);
		Real dis = 0;
		for (size_t i = 0; i < dim; ++i) {
			x[i] = rnd->normal.next();
			dis += x[i] * x[i];
		}
		dis = sqrt(dis);
		const auto r = rnd->uniform.nextNonStd(0.0, rcloud);
		for (size_t i = 0; i < dim; ++i) {
			variable()[i] = (best.variable()[i]) + r * x[i] / dis;
		}
		x.clear();
		CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
		return evaluate(pro, alg);
	}

	int IndDynDE::entropy(Real sigma, Problem *pro, Algorithm *alg, Random *rnd) {
		for (auto& i : variable())
			i = i + rnd->normal.nextNonStd(0, sigma);
		CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
		return evaluate(pro, alg);
	}
}