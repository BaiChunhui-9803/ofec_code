/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F05",
	"identifier": "CEC2015_EOP_F05",
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


#ifndef OFEC_F5_SHIFTED_ROTATED_KATSUURA_H
#define OFEC_F5_SHIFTED_ROTATED_KATSUURA_H

#include "../../global/classical/katsuura.h"

namespace ofec {
	namespace cec2015 {
		class F5_shifted_rotated_katsuura final : public katsuura
		{
		public:
			F5_shifted_rotated_katsuura(const ParameterMap &v);
			F5_shifted_rotated_katsuura(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_EOP_F05 = CEC2015::F5_shifted_rotated_katsuura;
}
#endif // ! OFEC_F5_SHIFTED_ROTATED_KATSUURA_H
