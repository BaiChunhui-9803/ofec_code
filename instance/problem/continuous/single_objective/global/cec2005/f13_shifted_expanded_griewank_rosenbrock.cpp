#include "f13_shifted_expanded_griewank_rosenbrock.h"

namespace ofec {
	namespace cec2005 {
		void ShiftedExpandedGriewankRosenbrock::addInputParameters() {
			m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		}

		void ShiftedExpandedGriewankRosenbrock::initialize_(Environment *env) {
			Function::initialize_(env);
			resizeVariable(m_number_variables);
			setDomain(-3, 1);
			setBias(-130);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F13");  //data path
		}

		void ShiftedExpandedGriewankRosenbrock::updateOptima(Environment *env) {
			setOriginalGlobalOpt();
			setGlobalOpt(m_translation.data());
			m_objective_accuracy = 1.0e-2;
		}

		void ShiftedExpandedGriewankRosenbrock::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
			Real result = 0;
			for (size_t i = 0; i < m_number_variables; ++i) {
				Real result_f2 = 0;
				Real result_f8 = 0;
				Real x_front = x[i] + 1;
				Real x_back = x[(i + 1) % m_number_variables] + 1;
				result_f2 += 100 * pow((x_back - x_front * x_front), 2.0) + (x_front - 1) * (x_front - 1);
				result_f8 += result_f2 * result_f2 / 4000.0 - cos(result_f2 / sqrt((Real)(i + 1))) + 1;
				result += result_f8;
			}
			result += m_bias;
			obj[0] = result;
		}
	}
}