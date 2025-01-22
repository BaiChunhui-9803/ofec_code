/********* Begin Register Information **********
{
	"name": "Ring-PSO",
	"identifier": "RingPSO",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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
// Created: 21 September 2011
// Modified: 21 Nov 2014
// Modified: 2 Dec 2019 by WangJunChen

/***********************************************************************************************
* @article{li2009niching,
  volume = {14},
  number = {1},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {150--169},
  year = {2009},
  author = {Xiaodong Li},
  title = {Niching without niching parameters: Particle swarm optimization using a ring topology}
}
****************************************************************************************************/

#ifndef OFEC_RINGPSO_H
#define OFEC_RINGPSO_H

#include "../../../../template/classic/particle_swarm_optimization/swarm.h"
#include "../../../../template/classic/particle_swarm_optimization/particle.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class RingParticle : public Particle {
	public:
		RingParticle(size_t num_obj, size_t num_con, size_t size_var) : Particle(num_obj, num_con, size_var) {}
		//RingParticle(const Solution<> & rhs) : Particle(rhs) {}
		void nextVelocity(const Solution<>* lbest, Real w, Real c1, Real c2, Random *rnd) override;
	};

	class RingSwarm : public Swarm<RingParticle> {
	public:
		enum class Topology { R2, R3, LHC_R2, LHC_R3 };
		RingSwarm(size_t size_pop, Topology topology, Environment *env);
		void setNeighborhood(Random *rnd) override;
	protected:
		bool m_is_neighbor_set;
		Topology m_topology;
	};

	class RingPSO : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(RingPSO)
	protected:
		std::unique_ptr<RingSwarm> m_pop;
		size_t m_pop_size;
		Real m_weight, m_accelerator1, m_accelerator2;
		RingSwarm::Topology m_topology;
		
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
		void initSwarm(Environment *env);
	};
}
#endif // !OFEC_RINGPSO_H
