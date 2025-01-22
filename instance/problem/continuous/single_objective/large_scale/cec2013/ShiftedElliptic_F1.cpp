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

#include "ShiftedElliptic_F1.h"
 
namespace ofec {
	namespace CEC2013 {
		ShiftedElliptic_F1::ShiftedElliptic_F1(const ParameterMap &v) :
			ShiftedElliptic_F1((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		ShiftedElliptic_F1::ShiftedElliptic_F1(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		ShiftedElliptic_F1::~ShiftedElliptic_F1() {
			delete[] mp_Ovector;
			delete[] mp_anotherz;
		}

		void ShiftedElliptic_F1::initialize() {
			m_variable_monitor = true;
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			ID = 1;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();

			setOriginalGlobalOpt();
			setGlobalOpt(mp_Ovector);
			m_initialized = true;
		}

		int ShiftedElliptic_F1::evaluateObjective(Real *x, std::vector<Real> &obj) {
			size_t i;
			for (i = 0; i < m_number_variables;++i) {
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}
			obj[0] = elliptic(mp_anotherz, m_number_variables);
			return kNormalEval;
		}
	}
}