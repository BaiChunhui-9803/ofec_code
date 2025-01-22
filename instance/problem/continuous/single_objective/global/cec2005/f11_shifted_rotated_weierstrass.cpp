#include "f11_shifted_rotated_weierstrass.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedWeierstrass::initialize_(Environment *env) {
			Weierstrass::initialize_(env);
			setConditionNumber(5);
			setBias(90);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F11");  //data path
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F11");
		}

		void ShiftedRotatedWeierstrass::updateOptima(Environment *env) {
			setOriginalGlobalOpt();
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}