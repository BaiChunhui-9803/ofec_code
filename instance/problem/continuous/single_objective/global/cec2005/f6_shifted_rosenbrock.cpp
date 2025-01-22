#include "f6_shifted_rosenbrock.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRosenbrock::initialize_(Environment *env) {
			Rosenbrock::initialize_(env);
			setDomain(-100, 100);
			setBias(390);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F06");  //data path
		}

		void ShiftedRosenbrock::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}