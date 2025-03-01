/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F08",
	"identifier": "CEC2015_EOP_F08",
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


#ifndef OFEC_F8_SHIFTED_ROTATED_GRIEWANK_ROSENBROCK_H
#define OFEC_F8_SHIFTED_ROTATED_GRIEWANK_ROSENBROCK_H

#include "../../global/classical/griewank_rosenbrock.h"

namespace ofec {
	namespace cec2015 {
		class F8_shifted_rotated_griewank_rosenbrock : public griewank_rosenbrock
		{
		public:
			F8_shifted_rotated_griewank_rosenbrock(const ParameterMap &v);
			F8_shifted_rotated_griewank_rosenbrock(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
		private:
		};
	}
	using CEC2015_EOP_F08 = CEC2015::F8_shifted_rotated_griewank_rosenbrock;
}
#endif // ! OFEC_F8_SHIFTED_ROTATED_GRIEWANK_ROSENBROCK_H

