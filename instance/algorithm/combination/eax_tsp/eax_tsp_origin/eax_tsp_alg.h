/********* Begin Register Information **********
{
	"name": "EAX-TSP",
	"identifier": "EAX_TSP",
	"tags": ["single-objective", "travelling salesman problem" ]
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
*  Created by Diao Yiya on 2023 11 13
*
*-------------------------------------------------------------------------------
*

*
*************************************************************************/

#ifndef OFEC_EAX_TSP_H
#define OFEC_EAX_TSP_H

#include "../../../../../core/algorithm/algorithm.h"
#include "environment.h"

namespace ofec {

#define CAST_EAX_TSP(alg) dynamic_cast<EAX_TSP*>(alg)

	class EAX_TSP : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(EAX_TSP)
	public:
		~EAX_TSP() = default;
		void addInputParameters() {}
		bool terminating() override {
			return Algorithm::terminating() || m_env->terminationCondition();
		}

		double bestObjective() const{
			return m_env->getfBestValue();
		}
		
	protected:
		std::unique_ptr<eax_tsp::TEnvironment> m_env;

		void initialize_(Environment *env) override;
		void run_(Environment* env) override;
	};
}

#endif // !OFEC_CMAES_H

