/********* Begin Register Information **********
{
	"name": "PMODE",
	"identifier": "PMODE",
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
* Wei Z, Gao W, Li G, et al. 
* A Penalty-Based Differential Evolution for Multimodal Optimization[J]. 
* IEEE Transactions on Cybernetics, 2021.
*-------------------------------------------------------------------------------
* Created on 14th April, 2022 by Junchen Wang
*********************************************************************************/

#ifndef OFEC_PMODE_H
#define OFEC_PMODE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../global/jade/jade_pop.h"

namespace ofec {
	class PenaltyIndDE : virtual public IndividualDE {
	public:
		PenaltyIndDE(size_t num_objs, size_t num_cons, size_t num_vars);
		int select(Environment *env) override;
		Real penalizedFitness() const { return m_penalized_fitness; }
		Real& penalizedFitness() { return m_penalized_fitness; }
	private:
		Real m_penalized_fitness;
	};

	class PMODE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(PMODE)
	protected:
		size_t m_pop_size;
		Real m_phi, m_epsilon, m_tilde_R, m_freq, m_tau, m_R_g;
		std::vector<Real> original_obj;
		std::list<Solution<>> m_S;
		PopJADE<PenaltyIndDE> m_pop;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		void updatePenaltyRadius(size_t g, const std::vector<size_t> &flag, Environment *env);
		void updateEliteSet(const PenaltyIndDE &best_ind, Real currentbest, Environment *env);

	public:
		Real getPenalizedFitness(const Solution<> &s, Environment *env) const;
		const std::list<Solution<>>& archiveSet() const { return m_S; }
		Real penaltyRadius() const { return m_R_g; }
	};
}

#endif // !OFEC_PMODE_H
