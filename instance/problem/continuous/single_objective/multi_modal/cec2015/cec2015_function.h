/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* class function defines basic operations for function/numerical optimization problems
*
*********************************************************************************/

#ifndef OFEC_CEC2015_FUNCTION_H
#define OFEC_CEC2015_FUNCTION_H

#include "../../../../../../core/problem/continuous/function.h"

namespace ofec {
	namespace CEC2015 {
		class CEC2015_function :public Function {
		public:
			CEC2015_function() = default;
			//CEC2015_function(const std::string &name, size_t size_var, size_t size_obj = 1);
			virtual void initialize_(Environment *env) override;
		protected:
			bool load_optima(const std::string &path);
			void load_optima_(const std::string &path);
			void set_optima();

			bool loadTranslation(const std::string &path);
			void load_translation_(const std::string &path);

			void evaluate_optima();
			void rotate(Real *x);
		protected:

		};
	}
}
#endif // !OFEC_CEC2015_FUNCTION_H
