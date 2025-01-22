#include "G03.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G03::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t i = 0; i < m_num_vars; i++) {
				m_domain.setRange(0., 1., i);
			}
			m_domain_update = true;
			for (size_t i = 0; i < m_num_vars; i++) {
				m_initial_domain.setRange(0., 1., i);
			}
			m_num_cons = 1;
			m_constraint.resize(1);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}
		}

		void G03::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			size_t i, j;
			Real t1 = 1.0, t2 = 0.0;
			obj[0] = 0.0;
			for (i = 0; i < 10; ++i)
			{
				t1 *= x[i];
			}
			obj[0] = -pow(sqrt(10),10)*t1;

			// evaluate constraint value
			for (i = 0; i < 10; ++i)
			{
				t2 += x[i]*x[i];
			}
			con[0] = t2 - 1-0.00000001;
			if (con[0] <= 0) con[0] = 0;

		}
	}
}

