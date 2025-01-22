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

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	namespace cec2015 {
		class Function : virtual public Continuous, public MetricsMMOP {
		protected:
			Real m_scale, m_bias;
			std::vector<Real> m_translation;
			std::vector<std::vector<Real>> m_rotation;

			virtual void initialize_(Environment *env) override;

			bool loadOptima(const std::string &path);
			void loadOptima_(const std::string &path);

			bool loadTranslation(const std::string &path);
			void loadTranslation_(const std::string &path);

			bool loadRotation(const std::string &path);
			void loadRotation_(const std::string &path);

			void evaluateOptima();
			void translate(Real *x);
			void scale(Real *x);
			void rotate(Real *x);
		};
	}
}
#endif // !OFEC_CEC2015_FUNCTION_H
