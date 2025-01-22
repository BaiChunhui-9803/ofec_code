/********* Begin Register Information **********
{
	"name": "ALC-PSO",
	"identifier": "ALC_PSO",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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
*********************************************************************************/

/******************************************************************************
* W.-N. Chen, J. Zhang, Y. Lin, N. Chen, Z.-H. Zhan, H. S.-H. Chung, Y. Li, and Y.-H. Shi, 
* Particle swarm optimization with an aging leader and challengers,
* IEEE transactions on evolutionary computation, vol. 17, no. 2, pp. 241¨C258, 2012.
*-------------------------------------------------------------------------------
* Created on 7th April, 2022 by Junchen Wang
*********************************************************************************/

#ifndef OFEC_ALC_PSO_H
#define OFEC_ALC_PSO_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/particle_swarm_optimization/particle.h"
#include "../../../../template/classic/particle_swarm_optimization/swarm.h"

namespace ofec {
	class ALC_PSO :virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(ALC_PSO)
	protected:
		size_t m_swarm_size, m_initial_lifespan, m_T;
		Real m_weight, m_accelerator1, m_accelerator2, m_pro;
		Swarm<Particle> m_swarm;
		Solution<> m_leader, m_challenger;
		size_t m_age, m_lifespan;
		std::array<Real, 3> m_pre_obj; // [0]: gbest, [1]: collective pbest, [2]: leader
		std::vector<Particle> m_old;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	private:
		void adjustLifespan();
		void generateChallenger(Environment *env);
		void evaluateChallenger(Environment *env);
	};
}

#endif // !OFEC_ALC_PSO_H
