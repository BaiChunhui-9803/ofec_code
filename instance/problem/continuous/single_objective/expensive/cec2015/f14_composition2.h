/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F14",
	"identifier": "CEC2015_EOP_F14",
	"problem tags": [ "EOP", "SOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#ifndef OFEC_F14_EXPENSIVE_COMPOSITION2_H
#define OFEC_F14_EXPENSIVE_COMPOSITION2_H

#include "../../global/cec2015/composition.h"

namespace ofec {
	namespace cec2015 {
		class F14_expensive_composition2 final: public composition_2015
		{
		public:
			F14_expensive_composition2(const ParameterMap &v);
			F14_expensive_composition2(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();

		private:

		};
	}
	using CEC2015_EOP_F14 = CEC2015::F14_expensive_composition2;
}
#endif // !OFEC_F14_EXPENSIVE_COMPOSITION2_H

