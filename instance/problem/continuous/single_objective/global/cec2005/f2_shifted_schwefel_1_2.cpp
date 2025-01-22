#include "f2_shifted_schwefel_1_2.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedSchwefel_1_2::initialize_(Environment *env) {
			Schwefel_1_2::initialize_(env);
			setBias(-450);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F02");  //data path
		}

		void ShiftedSchwefel_1_2::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-6;
		}
	}
}