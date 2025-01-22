#include "f14_shifted_rotated_expanded_scaffer_f6.h"
#include <numeric>

namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedExpandedScafferF6::initialize_(Environment *env) {
			ScafferF6::initialize_(env);
			setConditionNumber(3);
			setBias(-300);			
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F14");  //data path
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F14");		
		}

		void ShiftedRotatedExpandedScafferF6::updateOptima(Environment *env) {
			setOriginalGlobalOpt();
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}