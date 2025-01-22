#include "G21.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G21::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(0., 1000., 0);
			m_domain.setRange(0., 40., 1);
			m_domain.setRange(0., 40., 2);
			m_domain.setRange(100., 300., 3);
			m_domain.setRange(6.3, 6.7, 4);
			m_domain.setRange(5.9, 6.4, 5);
			m_domain.setRange(4.5, 6.25, 6);

			m_initial_domain.setRange(0., 1000., 0);
			m_initial_domain.setRange(0., 40., 1);
			m_initial_domain.setRange(0., 40., 2);
			m_initial_domain.setRange(100., 300., 3);
			m_initial_domain.setRange(6.3, 6.7, 4);
			m_initial_domain.setRange(5.9, 6.4, 5);
			m_initial_domain.setRange(4.5, 6.25, 6);
	
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 6;
			m_constraint.resize(6);
			m_constraint[0] = Constraint::kInequality;
			for (size_t i = 1; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}
		}

		void G21::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = x[0];

			// evaluate constraint value
			con[0] = -x[0] + 35 * (pow(x[1],0.6)) + 35 * (pow(x[2],0.6));
			con[1] = -300 * x[2] + 7500 * x[4] - 7500 * x[5] - 25 * x[3] * x[4] + 25 * x[3] * x[5] + x[2] * x[3] - 0.0001;
			con[2] = 100 * x[1] + 155.365 * x[3] + 2500 * x[6] - x[1] * x[3] - 25 * x[3] * x[6] - 15536.5 - 0.0001;
			con[3] = -x[4] + log(-x[3] + 900) - 0.0001;
			con[4] = -x[5] + log(x[3] + 300) - 0.0001;
			con[5] = -x[6] + log(-2 * x[3] + 700) - 0.0001;

			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;


		}
	}
}

