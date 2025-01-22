#include "G24.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G24::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(0., 3., 0);
			m_domain.setRange(0., 4., 1);
			m_domain_update = true;
			m_initial_domain.setRange(0., 3., 0);
			m_initial_domain.setRange(0., 4., 1);
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

		}

		void G24::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = -x[0] - x[1];

			// evaluate constraint value
			con[0] = -2 * pow(x[0],4) + 8 * pow(x[0],3) - 8 * pow(x[0],2) + x[1] - 2;
			con[1] = -4 * pow(x[0],4) + 32 * pow(x[0],3) - 88 * pow(x[0],2) + 96 * x[0] + x[1] - 36;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;


		}
	}
}

