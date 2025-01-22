#include "F4_shifted_rotated_schwefel.h"

namespace ofec {
	namespace cec2015 {
		F4_shifted_rotated_schwefel::F4_shifted_rotated_schwefel(const ParameterMap &v) :
			F4_shifted_rotated_schwefel((v.at("problem name")), (v.at("number of variables")), 1) {

			
		}
		F4_shifted_rotated_schwefel::F4_shifted_rotated_schwefel(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			schwefel(name, size_var, size_obj) {

			
		}

		void F4_shifted_rotated_schwefel::initialize() {
			m_variable_monitor = true;
			setDomain(-500, 500);
			setInitialDomain(-500, 500);
			std::vector<Real> v(m_number_variables, 420.9687);
			setOriginalGlobalOpt(v.data());
			
			m_variable_accuracy = 1.0e-2;
			setBias(400);

			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt(m_translation.data());
			m_initialized = true;
		}
		int F4_shifted_rotated_schwefel::evaluateObjective(Real *x, std::vector<Real> &obj) {
			return schwefel::evaluateObjective(x, obj);
		}
	}
}