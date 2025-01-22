#include "f9_shifted_rastrigin.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRastrigin::initialize_(Environment *env) {
			Rastrigin::initialize_(env);
			setDomain(-5, 5);
			setBias(-330);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F09");  //data path
		}

		void ShiftedRastrigin::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}