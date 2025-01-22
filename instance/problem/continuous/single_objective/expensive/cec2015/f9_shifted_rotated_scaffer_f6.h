/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F09",
	"identifier": "CEC2015_EOP_F09",
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


#ifndef OFEC_F9_SHIFTED_ROTATED_SCAFFER_F6_H
#define OFEC_F9_SHIFTED_ROTATED_SCAFFER_F6_H

#include "../../global/classical/scaffer_F6.h"

namespace ofec {
	namespace cec2015 {
		class F9_shifted_rotated_scaffer_F6 final : public scaffer_F6
		{
		public:
			F9_shifted_rotated_scaffer_F6(const ParameterMap &v);
			F9_shifted_rotated_scaffer_F6(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_EOP_F09 = CEC2015::F9_shifted_rotated_scaffer_F6;
}
#endif // ! OFEC_F9_SHIFTED_ROTATED_SCAFFER_F6_H

