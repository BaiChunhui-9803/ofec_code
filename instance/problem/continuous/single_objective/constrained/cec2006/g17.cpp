#include "G17.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G17::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(0., 400., 0);
			m_domain.setRange(0., 1000., 1);
			for (size_t j = 2; j < 4; j++) {
				m_domain.setRange(340., 420., j);
			}
			m_domain.setRange(-1000., 1000., 4);
			m_domain.setRange(0., 0.5236, 5);

			m_initial_domain.setRange(0., 400., 0);
			m_initial_domain.setRange(0., 1000., 1);
			for (size_t j = 2; j < 4; j++) {
				m_initial_domain.setRange(340., 420., j);
			}
			m_initial_domain.setRange(-1000., 1000., 4);
			m_initial_domain.setRange(0., 0.5236, 5);
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 4;
			m_constraint.resize(4);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}
		}

		void G17::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			Real y1 = 0.0, y2 = 0.0, t5 = 0.0;
			if (x[0] < 300)  y1 = 30 * x[0];
			else y1 = 31 * x[0];

			if (x[1] < 100)  y2 = 28 * x[1];
			else if (x[1] >= 100 && x[1] < 200)  y2 = 29 * x[1];
			else y2 = 30 * x[1];

			obj[0] = y1 + y2;

			// evaluate constraint value
			
			con[0] = -x[0] + 300 - (x[2] * x[3] / 131.078) * cos(1.48477 - x[5]) + 0.90798 * x[2] * x[2] * cos(1.47588) / 131.078;
			con[1] = -x[1] - x[2] * x[3] * cos(1.48477 + x[5]) / 131.078 + 0.90798 * x[3] * x[3] * cos(1.47588) / 131.078;
			con[2] = -x[4] - x[2] * x[3] * sin(1.48477 + x[5]) / 131.078 + 0.90798 * x[3] * x[3] * sin(1.47588) / 131.078;
			con[3] = 200 - x[2] * x[3] * sin(1.48477 - x[5]) / 131.078 + 0.90798 * x[2] * x[2] * sin(1.47588) / 131.078;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;


		}
	}
}

