#include "f10_shifted_rotated_rastrigin.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedRastrigin::initialize_(Environment *env) {
			Rastrigin::initialize_(env);
			setDomain(-5, 5);
			setBias(-330);
			setConditionNumber(2);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F10");  //data path
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F10");	
		}

		void ShiftedRotatedRastrigin::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}