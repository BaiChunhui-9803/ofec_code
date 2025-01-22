#include "G18.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G18::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < (m_num_vars-1); j++) {
				m_domain.setRange(-10., 10., j);
			}
			m_domain.setRange(0., 20., 8);
			m_domain_update = true;
			for (size_t j = 0; j < (m_num_vars - 1); j++) {
				m_initial_domain.setRange(-10., 10., j);
			}
			m_initial_domain.setRange(0., 20., 8);
			
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 13;
			m_constraint.resize(13);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
		}

		void G18::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {
			obj[0] = -0.5 * (x[0] * x[3] - x[1] * x[2] + x[2] * x[8] - x[4] * x[8] + x[4] * x[7] - x[5] * x[6]);

			// evaluate constraint value
			con[0] = x[2] * x[2] + x[3] * x[3] - 1;
			con[1] = x[8] * x[8] - 1;
			con[2] = x[4] * x[4] + x[5] * x[5] - 1;
			con[3] = x[0] * x[0] + (x[1] - x[8]) * (x[1] - x[8]) - 1;
			con[4] = (x[0] - x[4]) * (x[0] - x[4]) + (x[1] - x[5]) * (x[1] - x[5]) - 1;
			con[5] = (x[0] - x[6]) * (x[0] - x[6]) + (x[1] - x[7]) * (x[1] - x[7]) - 1;
			con[6] = (x[2] - x[4]) * (x[2] - x[4]) + (x[3] - x[5]) * (x[3] - x[5]) - 1;
			con[7] = (x[2] - x[6]) * (x[2] - x[6]) + (x[3] - x[7]) * (x[3] - x[7]) - 1;
			con[8] = x[6] * x[6] + (x[7] - x[8]) * (x[7] - x[8]) - 1;
			con[9] = x[1] * x[2] - x[0] * x[3];
			con[10] = -x[2] * x[8];
			con[11] = x[4] * x[8];
			con[12] = x[5] * x[6] - x[4] * x[7];

			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;
			if (con[6] <= 0) con[6] = 0;
			if (con[7] <= 0) con[7] = 0;
			if (con[8] <= 0) con[8] = 0;
			if (con[9] <= 0) con[9] = 0;
			if (con[10] <= 0) con[10] = 0;
			if (con[11] <= 0) con[11] = 0;
			if (con[12] <= 0) con[12] = 0;

		}
	}
}

