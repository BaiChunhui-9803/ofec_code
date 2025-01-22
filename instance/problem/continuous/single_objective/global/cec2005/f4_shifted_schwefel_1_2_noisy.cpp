#include "f4_shifted_schwefel_1_2_noisy.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedSchwefel_1_2_Noisy::initialize_(Environment *env) {
			Schwefel_1_2::initialize_(env);
			setBias(-450);	
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F04");  //data path		
		}

		void ShiftedSchwefel_1_2_Noisy::updateOptima(Environment *env) {
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-6;
		}

		void ShiftedSchwefel_1_2_Noisy::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
			Schwefel_1_2::evaluateOriginalObj(x, obj);
			obj[0] = (obj[0] - m_bias)*fabs(m_random->normal.next()) + m_bias;
		}
	}
}