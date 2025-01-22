/********* Begin Register Information **********
{
	"name": "MMOP_CEC2015_F14",
	"identifier": "CEC2015_MMOP_F14",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Problem Definitions and Evaluation Criteria for the CEC 2015
*  Competition on Single Objective Multi-Niche Optimization.
*******************************************************************************************/


#ifndef OFEC_CEC2015_F14_COMPOSITION6_H
#define OFEC_CEC2015_F14_COMPOSITION6_H

#include "../../expensive/cec2015/composition_2015.h"


namespace ofec {
	namespace cec2015 {
		class F14_composition2015_C6 final : public composition_2015
		{
		public:
			F14_composition2015_C6(const ParameterMap &v);
			F14_composition2015_C6(const std::string &name, size_t size_var, size_t size_obj);
			Real pre_opt_distance();
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();

			void rotate(size_t num, Real *x);
			void scale(size_t num, Real *x);

			bool loadTranslation(const std::string &path);
			void setTranslation();
			void set_weight(std::vector<Real>& weight, const std::vector<Real>&x);
		private:
			Real m_pre_opt_distance;
		};
	}
	using CEC2015_MMOP_F14 = cec2015::F14_composition2015_C6;
}
#endif // !OFEC_CEC2015_F14_COMPOSITION6_H










