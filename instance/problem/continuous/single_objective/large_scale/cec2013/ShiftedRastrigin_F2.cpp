/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
#include "ShiftedRastrigin_F2.h"
namespace ofec {
	namespace CEC2013 {
		ShiftedRastrigin_F2::ShiftedRastrigin_F2(const ParameterMap &v) :
			ShiftedRastrigin_F2((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		ShiftedRastrigin_F2::ShiftedRastrigin_F2(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		ShiftedRastrigin_F2::~ShiftedRastrigin_F2() {
			delete[] mp_Ovector;
			delete[] mp_anotherz;
		}

		void ShiftedRastrigin_F2::initialize() {
			m_variable_monitor = true;
			setDomain(-5, 5);
			setInitialDomain(-5, 5);
			ID = 2;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();

			setOriginalGlobalOpt();
			setGlobalOpt(mp_Ovector);
			m_initialized = true;
		}

		int ShiftedRastrigin_F2::evaluateObjective(Real *x, std::vector<Real> &obj) {
			size_t i;
			for (i = 0; i < m_number_variables; ++i) {
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}
			obj[0] = rastrigin(mp_anotherz, m_number_variables);
			return kNormalEval;
		}
	}
}