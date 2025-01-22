#include "G23.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G23::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(0., 300., 0);
			m_domain.setRange(0., 300., 1);
			m_domain.setRange(0., 100., 2);
			m_domain.setRange(0., 200., 3);
			m_domain.setRange(0., 100., 4);
			m_domain.setRange(0., 300., 5);
			m_domain.setRange(0., 100., 6);
			m_domain.setRange(0., 200., 7);
			m_domain.setRange(0.01, 0.03, 8);
			m_domain_update = true;

			m_initial_domain.setRange(0., 300., 0);
			m_initial_domain.setRange(0., 300., 1);
			m_initial_domain.setRange(0., 100., 2);
			m_initial_domain.setRange(0., 200., 3);
			m_initial_domain.setRange(0., 100., 4);
			m_initial_domain.setRange(0., 300., 5);
			m_initial_domain.setRange(0., 100., 6);
			m_initial_domain.setRange(0., 200., 7);
			m_initial_domain.setRange(0.01, 0.03, 8);
			
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 6;
			m_constraint.resize(6);
			m_constraint[0] = Constraint::kInequality;
			m_constraint[1] = Constraint::kInequality;
			for (size_t i = 2; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}
		}

		void G23::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = -9 * x[4] - 15 * x[7] + 6 * x[0] + 16 * x[1] + 10 * (x[5] + x[6]);

			// evaluate constraint value
			con[0] = x[8] * x[2] + 0.02 * x[5] - 0.025 * x[4];
			con[1] = x[8] * x[3] + 0.02 * x[6] - 0.015 * x[7];
			con[2] = x[0] + x[1] - x[2] - x[3];
			con[3] = 0.03 * x[0] + 0.01 * x[1] - x[8] * (x[2] + x[3]);
			con[4] = x[2] + x[5] - x[4];
			con[5] = x[3] + x[6] - x[7];
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;


		}
	}
}

