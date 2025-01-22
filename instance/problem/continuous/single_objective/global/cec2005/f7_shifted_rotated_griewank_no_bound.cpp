#include "f7_shifted_rotated_griewank_no_bound.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedGriewankNoBound::initialize_(Environment *env) {
			Griewank::initialize_(env);
			for (size_t i = 0; i < m_number_variables; i++)
				m_domain[i].limited = false;
			setDomainInitPop(0, 600.);
			setBias(-180);
			setConditionNumber(3);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F07");  //data path
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F07");
		}

		void ShiftedRotatedGriewankNoBound::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}