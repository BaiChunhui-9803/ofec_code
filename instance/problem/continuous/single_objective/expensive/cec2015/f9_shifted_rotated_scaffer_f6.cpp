#include "F9_shifted_rotated_scaffer_F6.h"

namespace ofec {
	namespace cec2015 {
		F9_shifted_rotated_scaffer_F6::F9_shifted_rotated_scaffer_F6(const ParameterMap &v) :
			F9_shifted_rotated_scaffer_F6((v.at("problem name")), (v.at("number of variables")), 1) {

			
		}
		F9_shifted_rotated_scaffer_F6::F9_shifted_rotated_scaffer_F6(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			scaffer_F6(name, size_var, size_obj) {

			
		}

		void F9_shifted_rotated_scaffer_F6::initialize() {
			m_variable_monitor = true;
			setDomain(-100, 100);
			setInitialDomain(-100, 100);
			setOriginalGlobalOpt();

			m_variable_accuracy = 1.0e-2;
			setBias(900);

			
			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt(m_translation.data());
			m_initialized = true;
		}
		int F9_shifted_rotated_scaffer_F6::evaluateObjective(Real *x, std::vector<Real> &obj) {
			return scaffer_F6::evaluateObjective(x, obj);
		}
	}
}