#include "G08.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G08::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(0., 10., j);
			}
			m_domain_update = true;
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(0., 10., j);
			}
			m_num_cons = 2;
			m_constraint.resize(2);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
		}

		void G08::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			Real a = pow((sin(2 * 3.14 * x[0])), 3);
			Real b = sin(2 * OFEC_PI * x[1]);

			if ( x[0] == 0 ) {
				x[0] = 1e-12;
			}
			obj[0] = -a * b / ((pow(x[0] ,3) * (x[0] + x[1])));
			
			con[0] = x[0] * x[0] - x[1] + 1;
			con[1] = 1 - x[0] + (x[1] - 4) * (x[1] - 4);
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;

		}
	}
}

