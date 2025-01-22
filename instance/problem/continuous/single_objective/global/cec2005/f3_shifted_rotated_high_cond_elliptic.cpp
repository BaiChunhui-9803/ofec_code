#include "f3_shifted_rotated_high_cond_elliptic.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedHighCondElliptic::initialize_(Environment *env) {
			Elliptic::initialize_(env);
			setBias(-450);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F03");  //data path
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F03");	
		}

		void ShiftedRotatedHighCondElliptic::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-6;
		}
	}
}