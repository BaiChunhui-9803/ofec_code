#include "F1_rotated_bent_cigar.h"

namespace ofec {
	namespace cec2015 {
		F1_rotated_bent_cigar::F1_rotated_bent_cigar(const ParameterMap &v) : F1_rotated_bent_cigar(v.at("problem name"), v.at("number of variables"), 1){
					
		}
		F1_rotated_bent_cigar::F1_rotated_bent_cigar(const std::string &name, size_t size_var, size_t size_obj ) :\
			problem(name, size_var, size_obj), bent_cigar(name, size_var, size_obj) {
			
		}

		void F1_rotated_bent_cigar::initialize() {
			m_variable_monitor = true;
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			
			setBias(100);
			setOriginalGlobalOpt();
			loadRotation("/instance/problem/continuous/expensive/cec2015");
			
			setGlobalOpt();
			m_initialized = true;
		}
		int F1_rotated_bent_cigar::evaluateObjective(Real *x, std::vector<Real> &obj) {
			bent_cigar::evaluateObjective(x, obj);
			return kNormalEval;
		}
	}
}