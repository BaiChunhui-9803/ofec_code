/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Yiya Diao & Junchen Wang & Changhe Li
* Email: diaoyiyacug@gmail.com & wangjunchen.chris@gmail.com & changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/

#ifndef OFEC_INSTANCE_H
#define OFEC_INSTANCE_H

#include "parameter/parameter_map.h"
#include "parameter/input_parameter_map.h"
#include "random/newran.h"

namespace ofec {
	class Instance  {
	public:
		virtual ~Instance() = default;

		/* Ready-only methods */
		const std::string& name() const { return m_name; }
		const std::shared_ptr<Random>& random() { return m_random; }
		const InputParameterMap& inputParameters() const { return m_input_parameters; }
		const ParameterMap& archivedParameters() const { return m_archived_parameters; }

		/* Write methods */
		void setName(const std::string &name) { m_name = name; }
		void setRandom(const std::shared_ptr<Random>& random) { m_random = random; }

		/**
		 * @brief_bch Records the input parameters into the archived parameters.
		 *
		 * This function records the current input parameters of the instance into the archived parameters.
		 * It first checks the validity of the input parameters and then outputs them to the archived parameters.
		 */
		void recordInputParameters();

		/**
		 * @brief_bch Restores the input parameters from the archived parameters.
		 *
		 * This function restores the input parameters of the instance from the archived parameters.
		 * It reads the parameters from the archived parameters and restores them to the input parameters.
		 */
		void restoreInputParameters();

	protected:
		Instance() = default;

		std::shared_ptr<Random> m_random;
		InputParameterMap m_input_parameters;
		ParameterMap m_archived_parameters;

	private:
		std::string m_name;
	};

#define OFEC_ABSTRACT_INSTANCE(identifier)\
protected:\
	 identifier() {\
		addInputParameters();\
	};\
private:

#define OFEC_CONCRETE_INSTANCE(identifier)\
public:\
	static identifier* create() {\
		return new identifier();\
	};\
OFEC_ABSTRACT_INSTANCE(identifier)
}

#endif