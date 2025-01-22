/********* Begin Register Information **********
{
	"name": "EAX-TSP2",
	"identifier": "EAX_TSP2",
	"problem tags": [ "ComOP", "SOP", "TSP" ]
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

#ifndef OFEC_EAX_TSP2_H
#define OFEC_EAX_TSP2_H

#include "../../../../../core/algorithm/algorithm.h"
#include "environment.h"

namespace ofec {
#define CAST_EAX_TSP2(alg) dynamic_cast<EAX_TSP2*>(alg)

	class EAX_TSP2 : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(EAX_TSP2)

	public:

		void setInputParameters() {}

		~EAX_TSP2() = default;

		bool terminating() override {
			return Algorithm::terminating() || m_env->terminationCondition();
		}

		eax_tsp2::TEnvironment& environment() {
			return *m_env.get();
		}

		double randomNumber() {
			return m_random->uniform.next();
		}

		Random* getRandom() {
			return m_random.get();
		}

	protected:
		std::unique_ptr<eax_tsp2::TEnvironment> m_env;

		void initialize_(Environment *env) override;
		void run_(Environment* env) override;
	};
}

#endif // !OFEC_EAX_TSP2_H

