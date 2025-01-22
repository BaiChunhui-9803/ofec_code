#include "G06.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {
		void G06::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(13., 100., 0);
			m_domain.setRange(0., 100., 1);
			m_domain_update = true;
			m_initial_domain.setRange(13., 100., 0);
			m_initial_domain.setRange(0., 100., 1);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
		}

		void G06::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			if (x[0] == 10) {
				x[0]=10+1e-12;
			}
			if (x[1] == 20) {
				x[1] = 20 + 1e-12;
			}
			obj[0] = pow((x[0]-10),3)+pow((x[1]-20),3);

			con[0] = -pow((x[0] - 5), 2) - pow((x[1] - 5), 2) + 100;
			con[1] = pow((x[0] - 6), 2) + pow((x[1] - 5), 2) - 82.81;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
		}
	}
}

