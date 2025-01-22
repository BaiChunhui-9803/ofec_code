/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F10",
	"identifier": "CEC2015_EOP_F10",
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


#ifndef OFEC_F10_HYBRID1_H
#define OFEC_F10_HYBRID1_H

#include "../../global/cec2015/hybrid.h"

namespace ofec {
	namespace cec2015 {
		class F10_hybrid1 final : public hybrid
		{
		public:
			F10_hybrid1(const ParameterMap &v);
			F10_hybrid1(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();
		private:
		
		};
	}
	using CEC2015_EOP_F10 = CEC2015::F10_hybrid1;
}
#endif // ! OFEC_F10_HYBRID1_H


