#include "F6_shifted_rotated_happy_cat.h"

namespace ofec {
	namespace cec2015 {
		F6_shifted_rotated_happy_cat::F6_shifted_rotated_happy_cat(const ParameterMap &v) :
			F6_shifted_rotated_happy_cat((v.at("problem name")), (v.at("number of variables")), 1) {

			
		}
		F6_shifted_rotated_happy_cat::F6_shifted_rotated_happy_cat(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			happy_cat(name, size_var, size_obj) {

			
		}

		void F6_shifted_rotated_happy_cat::initialize() {
			m_variable_monitor = true;
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			setOriginalGlobalOpt();

			setBias(600);

			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt(m_translation.data());
			m_initialized = true;
		}
		int F6_shifted_rotated_happy_cat::evaluateObjective(Real *x, std::vector<Real> &obj) {
			return happy_cat::evaluateObjective(x, obj);
		}
	}
}