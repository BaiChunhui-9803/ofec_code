#include "ATF3.h"
#include <math.h>


namespace ofec {
	namespace ATF {
		void _3::addInputParameters() {
			m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		}
		void _3::initialize_(Environment* env) {
			Function::initialize_(env);
			//m_variable_monitor = true;
			//auto& v = GET_PARAM(m_id_param);
			//resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(-10., 5., 0);
			m_domain.setRange(-5., 5., 1);
			//m_domain_update = true;
			//m_initial_domain.setRange(-10., 5., 0);
			//m_initial_domain.setRange(-5., 5., 1);
			m_number_constraints = 1;
			m_constraint.resize(1);
			m_constraint[0] = Constraint::kInequality;
		}

		void _3::evaluateObjectiveAndConstraint(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = x[0] * x[0] - 10 * cos(2 * 3.14 * x[0]) + x[1] * x[1] - 10 * cos(2 * 3.14 * x[1]) + 20;

			con[0] = 3 * (x[0] + 9) * (x[0] + 9) + x[1] * x[1] - 2;

			if (con[0] <= 0) con[0] = 0;

		}
	}
}