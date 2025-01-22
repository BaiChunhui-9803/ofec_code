/********* Begin Register Information **********
{
	"name": "SDCMOP",
	"identifier": "SDCMOP",
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

// Created: 9 June 2020 by Qingshan Tan


#ifndef SDCMOP_H
#define SDCMOP_H

#include "../DMOPs.h"
//#include "../../dynamic/dynamic.h"

namespace ofec {
	class SDCMOP :public DMOPs {
	public:
		SDCMOP(const ParameterMap &v);
		SDCMOP(const std::string &name, size_t size_var, size_t size_obj, size_t type, size_t size_peak);
		~SDCMOP() {};
		enum class dynamic_factor {
			Normal,
			Change_pattern, Change_PF_curvature, Change_PF_seg, Change_num_objective, Change_num_pri, Change_num_pub,
			Change_location_peak, Change_width_peak, Change_num_PS, Change_center_PS, Change_radius_PS, Change_detect, Change_predict
		};
		void initialize();
		void generateAdLoadPF();
		void update_problem(Real t);
		void set_change_type(size_t i);
		dynamic_factor get_dynamic_type() { return dynamic_type; }
		/*void set_updated_state(bool it) { updated_state = it; }
		bool get_updated_state() { return updated_state; }*/
		void get_1d_sample_PF(size_t n);
	protected:
		void evaluate_obj_nd_con(Real *x, std::vector<Real>& obj, std::vector<Real> &con);
		//bool updated_state = false;//denote if updated problem
		dynamic_factor dynamic_type;
	private:
		Real m_T = 1. / 4;//denote the period of the PF
		Real m_alpha = 1. / 5;//denote the curvature of the PF
		size_t m_num_PS;//the number of PS
		size_t m_num_pri = 4;//the number of private variables
		size_t m_num_pub = 4;//the number of public variables
		std::vector<Real> m_R;//the radius of the PS
		std::vector<Real> m_R0;//the initial radius of the PS
		std::vector<std::vector<Real>> m_H;//store the peaks of all subobjectives
		std::vector<std::vector<Real>> m_H0;//store the initial value of H
		std::vector<std::vector<Real>> m_W;//store the width of peaks of all subobjectives
		std::vector<std::vector<Real>> m_W0;//store the initial value of W
		std::vector<std::vector<Real>> m_center_of_PS;//store the location of the center of PS
		std::vector<std::vector<Real>> m_center_of_PS0;//store the initial location of the center of PS
		std::vector<std::vector<std::vector<Real>>> m_variables;//store the location of peaks of all subobjectives
		std::vector<std::vector<std::vector<Real>>> m_x0;//store the initial value of m_variables
		std::vector<std::vector<Real>> m_sample_PF;//�������M-1����Ŀ��Ĳ�����
	};
}
#endif //SDMOP_H
