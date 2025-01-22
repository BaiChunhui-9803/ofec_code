#include "f6_sr_expanded_himmelblau.h"

namespace ofec {
	namespace cec2015 {
		//F6_SR_expanded_himmelblau::F6_SR_expanded_himmelblau(const ParameterMap &v) :
		//	F6_SR_expanded_himmelblau((v.at("problem name")), (v.at("number of variables")), 1) {
		//	
		//}
		//F6_SR_expanded_himmelblau::F6_SR_expanded_himmelblau(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
		//	CEC2015_function(name, size_var, size_obj) {
		//	
		//}

		void F6_SR_expanded_himmelblau::initialize_(Environment *env) {
			Function::initialize_(env);
			setDomain(-100, 100);
			setInitialDomain(-100, 100);

			m_variable_accuracy = 0.01;
			m_objective_accuracy = 1.e-4;



			m_bias = 600.;
			m_scale = 5;


			loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
			loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");

			// 2^Dim gopt
			evaluateOptima();
		}
		void F6_SR_expanded_himmelblau::evaluateObjective(Real *x, std::vector<Real> &obj) {

			size_t i;

			for (i = 0; i < m_number_variables; ++i)
				x[i] -= m_translation[i];
			scale(x);
			rotate(x);

			obj[0] = 0.0;

			for (i = 0; i<m_number_variables - 1; i = i + 2)
			{
				x[i] += 3.0;
				x[i + 1] += 2.0;//shift to orgin
				obj[0] += pow((x[i] * x[i] + x[i + 1] - 11.0), 2.0) + pow((x[i] + x[i + 1] * x[i + 1] - 7.0), 2.0);
			}
			obj[0] += m_bias;
			//return kNormalEval;
		}
	}
}