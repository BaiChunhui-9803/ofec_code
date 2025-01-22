/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F06",
	"identifier": "CEC2015_EOP_F06",
	"problem tags": [ "EOP", "SOP", "ConOP" ]
}
*********** End Register Information **********/

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


#ifndef OFEC_F6_SHIFTED_ROTATED_HAPPY_CAT_H
#define OFEC_F6_SHIFTED_ROTATED_HAPPY_CAT_H

#include "../../global/classical/happy_cat.h"

namespace ofec {
	namespace cec2015 {
		class F6_shifted_rotated_happy_cat final : public happy_cat
		{
		public:
			F6_shifted_rotated_happy_cat(const ParameterMap &v);
			F6_shifted_rotated_happy_cat(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_EOP_F06 = CEC2015::F6_shifted_rotated_happy_cat;
}
#endif // ! OFEC_F6_SHIFTED_ROTATED_HAPPY_CAT_H

