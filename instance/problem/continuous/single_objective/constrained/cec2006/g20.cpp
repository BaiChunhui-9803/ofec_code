#include "G20.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G20::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(0., 10., j);
			}
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(0., 10., j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 20;
			m_constraint.resize(20);
			for (size_t i = 0; i < 6; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
			for (size_t i = 6; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}
		}

		void G20::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			size_t i, j;
			Real t3 = 0.0, t4 = 1.0, t5 = 0.0;
			obj[0] = 0.0;
			for (int i = 0; i < 20; ++i)
			{
				t3 += pow(cos(x[i]), 4);
			}
			for (int i = 0; i < 20; ++i)
			{
				t4 *= pow(cos(x[i]), 2);
			}
			for (int i = 0; i < 20; ++i)
			{
				t5 += (i + 1) * pow(x[i], 2);
			}
			obj[0] = -abs((t3 - 2 * t4) / sqrt(t5));

			// evaluate constraint value
			Real t1 = 1.0, t2 = 0.0;
			for (int i = 0; i < 20; ++i)
			{
				t1 *= x[i];
			}
			con[0] = 0.75 - t1;
			for (int i = 0; i < 20; ++i)
			{
				t2 += x[i];
			}
			con[1] = t2 - 7.5 * 20;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;


		}
	}
}

