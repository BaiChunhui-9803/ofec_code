/********* Begin Register Information **********
{
	"name": "DHNDE",
	"identifier": "DHNDE",
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
@article{wang2022multimodal,
  volume = {238},
  journal = {Knowledge-Based Systems},
  pages = {107972},
  year = {2022},
  author = {Kai Wang and Wenyin Gong and Libao Deng and Ling Wang},
  title = {Multimodal optimization via dynamically hybrid niching differential evolution}
}
*-------------------------------------------------------------------------------
* Created on 12th April, 2022 by Junchen Wang
*********************************************************************************/

#ifndef OFEC_DHNDE_H
#define OFEC_DHNDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class IndDHNDE : public IndividualDE {
	protected:
		size_t m_counter = 0;
	public:
		void resetCounter() { m_counter = 0; }
		void increaseCounter() { m_counter++; }
		size_t counter() const { return m_counter; }
		IndDHNDE(Environment *env) : IndividualDE(env) {}
		IndDHNDE(size_t number_objectives, size_t num_cons, size_t num_vars) : IndividualDE(number_objectives, num_cons, num_vars) {}
		IndDHNDE(const Solution<> &sol) : IndividualDE(sol) {}
	};

	class DHNDE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(DHNDE)
	protected:
		size_t m_pop_size, m_m, m_max_io, m_max_s, m_max_T;
		Real m_scaling_factor, m_Cr, m_d_min, m_f_min;
		de::MutateStrategy m_ms;
		PopulationDE<IndDHNDE> m_pop;
		std::list<Solution<>> m_A_io; // archive for inferior offspring
		std::list<Solution<>> m_A_os; // archive for optimal solutions
		bool m_without_crowding;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	private:
		void archiveUpdate(const Solution<> &x, Real d, Environment *env);
		void CDE_Aio(Environment *env);
		void INSDE(Environment *env);
	};
}

#endif // !OFEC_DHNDE_H
