/*************************************************************************
* Project: Library of Open Franmeworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Operator adaptator of genetic learning algorithm
*
*********************************************************************************/

#ifndef OFEC_GL_ADAPTER_H
#define OFEC_GL_ADAPTER_H
#include <vector>
#include <memory>
#include "../../../../../core/definition.h"
#include "../../../../../core/problem/problem.h"
#include "../../../../../core/algorithm/algorithm.h"
#include <list>

namespace ofec {
	template<typename TInd>
	class PopGL;

	template<typename TInd>
	class AdaptorGL {
	public:
		AdaptorGL(Real alpha) : m_alpha(alpha) {}
		AdaptorGL() = default;
		virtual ~AdaptorGL() {}
		virtual void updateProbability(Environment *env,
			PopGL<TInd> &pop,
			const std::vector<Real> &weight,
			const std::vector<int> *index = nullptr) {}
		virtual void updateProbability(Environment *env,
			std::vector<TInd> &pop,
			const std::vector<Real> &weight,
			const std::vector<int> *index = nullptr) {};
		virtual void createSolution(Environment *env, Random *rnd,
			PopGL<TInd> &pop,
			std::vector<TInd> &offspring) = 0;
		virtual int updateSolution(Environment *env,
			PopGL<TInd> &pop,
			std::vector<TInd> &offspring, int &num_improve) = 0;
		virtual void printInfo() {}
		Real getAlpha() const { return m_alpha; }
		void setAlpha(Real alpha) { m_alpha = alpha; }
	protected:
		Real m_alpha = 0.8;
	};
}

#endif //OFEC_GL_ADAPTER_H