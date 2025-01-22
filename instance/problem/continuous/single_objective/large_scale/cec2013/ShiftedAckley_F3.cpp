/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 8 August 2017
// Last modified:
#include "ShiftedAckley_F3.h"
namespace ofec {
	namespace CEC2013 {
		ShiftedAckley_F3::ShiftedAckley_F3(const ParameterMap &v) : 
			ShiftedAckley_F3((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		ShiftedAckley_F3::ShiftedAckley_F3(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		ShiftedAckley_F3::~ShiftedAckley_F3() {
			delete[] mp_Ovector;
			delete[] mp_anotherz;
		}

		void ShiftedAckley_F3::initialize() {
			m_variable_monitor = true;
			setDomain(-32, 32);
			setInitialDomain(-32, 32);
			ID = 3;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();

			setOriginalGlobalOpt();
			setGlobalOpt(mp_Ovector);
			m_initialized = true;
		}

		int ShiftedAckley_F3::evaluateObjective(Real *x, std::vector<Real> &obj) {
			size_t    i;
			for (i = 0; i < m_number_variables; ++i) {
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}
			obj[0] = ackley(mp_anotherz, m_number_variables);
			return kNormalEval;
		}
	}
}