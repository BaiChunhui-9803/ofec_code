#include "f4_sr_expanded_decreasing_minima.h"

namespace ofec {
	namespace cec2015 {
		//F4_SR_expanded_decreasing_minima::F4_SR_expanded_decreasing_minima(const ParameterMap &v) :
		//	F4_SR_expanded_decreasing_minima((v.at("problem name")), (v.at("number of variables")), 1) {
		//	
		//}
		//F4_SR_expanded_decreasing_minima::F4_SR_expanded_decreasing_minima(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
		//	CEC2015_function(name, size_var, size_obj) {
		//	
		//}

		void F4_SR_expanded_decreasing_minima::initialize_(Environment *env) {
			Function::initialize_(env);
			setDomain(-100, 100);
			setInitialDomain(-100, 100);
			
			m_variable_accuracy = 0.01;
			m_objective_accuracy = 1.e-4;



			m_bias = 400.;
			m_scale = 20;


			loadOptima("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
			loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2015/");

			// 5Dim  : 1 gopt and 15 lopt 
			// 10Dim : 1 gopt and 55 lopt 
			// 20Dim : 1 gopt and 210 lopt 
			evaluateOptima();



		}
		void F4_SR_expanded_decreasing_minima::evaluateObjective(Real *x, std::vector<Real> &obj) {

			size_t i;

			for (i = 0; i < m_number_variables; ++i)
				x[i] -= m_translation[i];
			scale(x);
			rotate(x);

			obj[0] = 0.0;
			
			for (i = 0; i<m_number_variables; i++)
			{
				x[i] += 0.1;
				if ((x[i] <= 1.0)&(x[i] >= 0.0))
				{
					obj[0] += -exp(-2.0*log(2.0)*pow((x[i] - 0.1) / 0.8, 2.0))*pow(sin(5.0*OFEC_PI*x[i]), 6.0);
				}
				else
				{
					obj[0] += pow(x[i], 2.0);
				}
			}
		
			obj[0] += 1.0*m_number_variables;
			obj[0] += m_bias;
	//	return kNormalEval;
		}
	}
}