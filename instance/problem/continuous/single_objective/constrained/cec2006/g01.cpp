#include "G01.h"
#include "../../../../../core/instance_manager.h"


namespace ofec {
	namespace CEC2006 {

		void G01::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < 9; j++) {
				m_domain.setRange(0., 1., j);
			}
			for (size_t j = 9; j < 12; j++) {
				m_domain.setRange(0., 100., j);
			}
			for (size_t j = 12; j < 13; j++) {
				m_domain.setRange(0., 1., j);
			}
			m_domain_update = true;
			for (size_t j = 0; j < 9; j++) {
				m_initial_domain.setRange(0., 1., j);
			}
			for (size_t j = 9; j < 12; j++) {
				m_initial_domain.setRange(0., 100., j);
			}
			for (size_t j = 12; j < 13; j++) {
				m_initial_domain.setRange(0., 1., j);
			}
			m_num_cons = 9;
			m_constraint.resize(9);
			for (size_t i = 0;i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
		}

		void G01::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			size_t i, j;
			Real t1 = 0.0, t2 = 0.0, t3 = 0.0;
			obj[0] = 0.0;
			for (i = 0; i < 4; ++i)
			{
				t1 += x[i];
			}
			for (i = 0; i < 4; ++i)
			{
				t2 += x[i]*x[i];
			}
			for (i = 4; i < 13; ++i)
			{
				t3 += x[i];
			}
			obj[0] = 5*t1-5*t2-t3;

			// evaluate constraint value

			con[0] = 2 * x[0] + 2 * x[1] + x[9] + x[10] - 10;
			con[1] = 2 * x[0] + 2 * x[2] + x[9] + x[11] - 10;
			con[2] = 2 * x[1] + 2 * x[2] + x[10] + x[11] - 10;
			con[3] = -8 * x[0] + x[9];
			con[4] = -8 * x[1] + x[10];
			con[5] = -8 * x[2] + x[11];
			con[6] = -2 * x[3] - x[4] + x[9];
			con[7] = -2 * x[5] - x[6] + x[10];
			con[8] = -2 * x[7] - x[8] + x[11];
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;
			if (con[6] <= 0) con[6] = 0;
			if (con[7] <= 0) con[7] = 0;
			if (con[8] <= 0) con[8] = 0;

		}
	}
}

