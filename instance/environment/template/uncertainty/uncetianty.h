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
*
*********************************************************************************/

#ifndef OFEC_UNCERTIANTY_H
#define OFEC_UNCERTIANTY_H

#include"dynamic.h"
#include"noisy.h"

namespace ofec {
#define GET_NOISY(pro) dynamic_cast<Uncertianty*>(pro)

	class Uncertianty : virtual public Noisy,virtual public Dynamic {
	private:
		virtual void getDynamicObjective(const SolutionBase& s, std::vector<Real>& objs, int t = 0) override {}
		virtual void getNoisyObjective(const SolutionBase& s, std::vector<Real>& objs, Random *rnd = nullptr) override {}

	public:
		virtual void getNoisyObjective(const SolutionBase& s, std::vector<Real>& objs, Random *rnd= nullptr) override {}
		// to do
		virtual int updateEvaluationTag(SolutionBase& s, Algorithm* alg) override {}

	protected:
		// to do
		virtual void initialize_()override {}



	};
}

#endif