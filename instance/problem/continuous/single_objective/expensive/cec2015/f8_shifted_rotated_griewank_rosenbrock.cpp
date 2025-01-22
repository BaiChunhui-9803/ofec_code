#include "F8_shifted_rotated_griewank_rosenbrock.h"

namespace ofec {
	namespace cec2015 {
		F8_shifted_rotated_griewank_rosenbrock::F8_shifted_rotated_griewank_rosenbrock(const ParameterMap &v) :
			F8_shifted_rotated_griewank_rosenbrock((v.at("problem name")), (v.at("number of variables")), 1) {

			
		}
		F8_shifted_rotated_griewank_rosenbrock::F8_shifted_rotated_griewank_rosenbrock(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			griewank_rosenbrock(name, size_var, size_obj) {

			
		}

		void F8_shifted_rotated_griewank_rosenbrock::initialize() {
			m_variable_monitor = true;
			setDomain(-5, 5);
			setInitialDomain(-5., 5.);
			setOriginalGlobalOpt();
			setBias(800);

			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt(m_translation.data());
			m_initialized = true;
		}
		int F8_shifted_rotated_griewank_rosenbrock::evaluateObjective(Real *x, std::vector<Real> &obj) {
			return griewank_rosenbrock::evaluateObjective(x, obj);
		}
	}
}