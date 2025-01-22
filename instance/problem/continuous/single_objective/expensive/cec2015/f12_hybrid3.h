/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F12",
	"identifier": "CEC2015_EOP_F12",
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


#ifndef OFEC_F12_HYBRID3_H
#define OFEC_F12_HYBRID3_H

#include "../../global/cec2015/hybrid.h"

namespace ofec {
	namespace cec2015 {
		class F12_hybrid3 final : public hybrid
		{
		public:
			F12_hybrid3(const ParameterMap &v);
			F12_hybrid3(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();
		private:
			
		};
	}
	using CEC2015_EOP_F12 = CEC2015::F12_hybrid3;
}
#endif // ! OFEC_F12_HYBRID3_H




