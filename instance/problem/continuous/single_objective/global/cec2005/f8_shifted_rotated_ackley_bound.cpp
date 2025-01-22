#include "f8_shifted_rotated_ackley_bound.h"


namespace ofec {
	namespace cec2005 {
		void ShiftedRotatedAckleyBound::initialize_(Environment *env) {
			Ackley::initialize_(env);
			setDomain(-32, 32);
			setBias(-140);
			setConditionNumber(100);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F08");  //data path
			for (size_t j = 1; j < m_number_variables / 2 + 1; ++j) {
				m_translation[(2 * j - 1) - 1] = -32;
				m_translation[2 * j - 1] = m_random->uniform.nextNonStd<Real>(-32, 32);
			}
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/");	
		}

		void ShiftedRotatedAckleyBound::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}
	}
}