#include "F3_shifted_rotated_weierstrass.h"

namespace ofec {
	namespace cec2015 {
		F3_shifted_rotated_weierstrass::F3_shifted_rotated_weierstrass(const ParameterMap &v) :F3_shifted_rotated_weierstrass((v.at("problem name")), (v.at("number of variables")), 1) {

			
		}
		F3_shifted_rotated_weierstrass::F3_shifted_rotated_weierstrass(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			weierstrass(name, size_var, size_obj) {

		
		}

		void F3_shifted_rotated_weierstrass::initialize() {
			setDomain(-0.5, 0.5);
			setInitialDomain(-0.5, 0.5);
			m_variable_monitor = true;
			setBias(300);
			setOriginalGlobalOpt();
			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt(m_translation.data());
			m_initialized = true;
		}
		int F3_shifted_rotated_weierstrass::evaluateObjective(Real *x, std::vector<Real> &obj) {
			return weierstrass::evaluateObjective(x, obj);
		}
	}
}