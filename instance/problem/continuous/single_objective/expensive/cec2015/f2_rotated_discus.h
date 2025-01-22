/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F02",
	"identifier": "CEC2015_EOP_F02",
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


#ifndef OFEC_F2_ROTATED_DISCUS_H
#define OFEC_F2_ROTATED_DISCUS_H

#include "../../global/classical/discus.h"

namespace ofec {
	namespace cec2015 {
		class F2_rotated_discus final : public discus
		{
		public:
			F2_rotated_discus(const ParameterMap &v);
			F2_rotated_discus(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_EOP_F02 = CEC2015::F2_rotated_discus;
}
#endif // ! OFEC_F2_ROTATED_DISCUS_H
