/********* Begin Register Information **********
{
	"name": "Multi-DOP",
	"identifier": "Multi_DOP",
	"problem tags": [ "DMOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Qingshan Tan
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

* see https://github.com/Changhe160/OFEC for more information
*************************************************************************/

/************************************************************************

************************************************************************/

// Created: 27 June 2020 by Qingshan Tan


#ifndef MULTI_DOP_H
#define MULTI_DOP_H

#include "../DMOPs.h"
//#include "../../dynamic/dynamic.h"

namespace ofec {
	class Multi_DOP :public DMOPs {
	public:
		Multi_DOP(const ParameterMap &v);
		Multi_DOP(const std::string &name, size_t size_var, size_t size_obj, size_t type, size_t size_peak);
		~Multi_DOP() {};
		enum class dynamic_factor {
			Normal,
			Change_pattern, Change_num_objective, Change_num_pri, Change_num_pub,
			Change_location_peak, Change_width_peak, Change_detect, Change_predict
		};
		void initialize();
		void generateAdLoadPF() {}
		void update_problem(Real t);
		void set_change_type(size_t i);
		dynamic_factor get_dynamic_type() { return dynamic_type; }
		void get_1d_sample_PF(size_t n) {}
	protected:
		void evaluateObjective(Real *x, std::vector<Real> &obj);
		dynamic_factor dynamic_type;
	private:
		size_t m_num_pri = 4;//the number of private variables
		size_t m_num_pub = 4;//the number of public variables
		std::vector<std::vector<Real>> m_H;//store the peaks of all subobjectives
		std::vector<std::vector<Real>> m_H0;//store the initial value of H
		std::vector<std::vector<Real>> m_W;//store the width of peaks of all subobjectives
		std::vector<std::vector<Real>> m_W0;//store the initial value of W
		std::vector<std::vector<std::vector<Real>>> m_variables;//store the location of peaks of all subobjectives
		std::vector<std::vector<std::vector<Real>>> m_x0;//store the initial value of m_variables
	};
}
#endif //Multi_DOP_H
