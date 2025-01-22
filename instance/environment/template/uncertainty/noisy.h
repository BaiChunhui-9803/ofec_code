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

// Created: 9 June 2021
// Last modified:

#ifndef OFEC_NOISY_H
#define OFEC_NOISY_H
#include "../problem.h"

namespace ofec {
#define GET_NOISY dynamic_cast<Noisy*>(pro)

	class Noisy : virtual public Problem {
	protected:
		// there are three sources of noisy
		bool m_flag_noisy_from_objective = false;
		Real m_objective_noisy_severity = 0.0;
		bool m_flag_noisy_from_variable = false;
		Real m_variable_noisy_severity = 0.0;
		bool m_flag_noisy_from_environment = false;
		Real m_environment_noisy_severity = 0.0;

	public:
		Noisy();
		virtual ~Noisy() = default;
		void copy(const Noisy &noisy);

		bool flagNoisy() const{
			return m_flag_noisy_from_objective|| m_flag_noisy_from_variable || m_flag_noisy_from_environment;
		}

		void setFlagNoisyFromObjective(bool flag_noisy_obj) {
			m_flag_noisy_from_objective = flag_noisy_obj;
		}

		bool getFlagNoisyFromObjective()const {
			return m_flag_noisy_from_objective;
		}

		void setFlagNoisyFromVariable(bool flag_noisy_var) {
			m_flag_noisy_from_variable = flag_noisy_var;
		}

		bool getFlagNoisyFromVariable()const {
			return m_flag_noisy_from_variable;
		}

		void setVariableNoisyServerity(Real val) {
			m_variable_noisy_severity = val;
		}

		void setFlagNoisyFromEnvironment(bool flag_noisy_env) {
			m_flag_noisy_from_environment = flag_noisy_env;
		}

		bool getFlagNoisyFromEnvironment()const {
			return m_flag_noisy_from_environment;
		}

		virtual void getNoisyObjective(const SolutionBase& s, std::vector<Real>& objs, Random *rnd) = 0;

	protected:
		void initialize_() override;
	};
}

#endif