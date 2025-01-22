#include "G19.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G19::initialize_() {
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
			m_num_cons = 5;
			m_constraint.resize(5);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
		}

		void G19::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			Real t1 = 0.0,t2 = 0.0,t3 = 0.0, t4 = 0.0, t5 = 0.0,t6 = 0.0,t7 = 0.0, t8 = 0.0,t9 = 0.0, t10 = 0.0,t11 = 0.0,t12=0.0;
			Real a[10][5] = { {-16,2,0,1,0},{0,-2,0,0.4,2},{-3.5,0,2,0,0},{0,-2,0,-4,-1},{0,-9,-2,1,-2.8},{2,0,-4,0,0},{-1,-1,-1,-1,-1},{-1,-2,-3,-2,-1},{1,2,3,4,5},{1,1,1,1,1} };
			Real b[10] = { -40,-2,-0.25,-4,-4,-1,-40,-60,5,1 };
			Real c[5][5] = { {30,-20,-10,32,-10},{-20,39,-6,-31,32},{-10,-6,10,-6,-10},{32,-31,-6,39,-20},{-10,32,-10,-20,30} };
			Real d[5] = { 4,8,10,6,2 };
			Real e[5] = { -15,-27,-36,-18,-12 };
			for (int i = 0; i < 5; ++i)
			{
				t1 += c[i][0] * x[10 + i];
			}
			for (int i = 0; i < 5; ++i)
			{
				t2 += c[i][1] * x[10 + i];
			}
			for (int i = 0; i < 5; ++i)
			{
				t3 += c[i][2] * x[10 + i];
			}
			for (int i = 0; i < 5; ++i)
			{
				t4 += c[i][3] * x[10 + i];
			}
			for (int i = 0; i < 5; ++i)
			{
				t5 += c[i][4] * x[10 + i];
			}
			for (int i = 0; i < 10; ++i)
			{
				t6 += a[i][0] * x[i];
			}
			for (int i = 0; i < 10; ++i)
			{
				t7 += a[i][1] * x[i];
			}
			for (int i = 0; i < 10; ++i)
			{
				t8 += a[i][2] * x[i];
			}
			for (int i = 0; i < 10; ++i)
			{
				t9 += a[i][3] * x[i];
			}
			for (int i = 0; i < 10; ++i)
			{
				t10 += a[i][4] * x[i];
			}
			for (int i = 0; i < 5; ++i)
			{
				t11 += d[i] * (x[10 + i] * x[10 + i] * x[10 + i]);
			}
			for (int i = 0; i < 10; ++i)
			{
				t12 += b[i] * x[i];
			}

			obj[0] = t1 * x[10] + t2 * x[11] + t3 * x[12] + t4 * x[13] + t5 * x[14] + 2 * t11 - t12;

			// evaluate constraint value
			con[0] = -2 * t1 - 3 * d[0] * (x[10] * x[10]) - e[0] + t6;
			con[1] = -2 * t2 - 3 * d[1] * (x[11] * x[11]) - e[1] + t7;
			con[2] = -2 * t3 - 3 * d[2] * (x[12] * x[12]) - e[2] + t8;
			con[3] = -2 * t4 - 3 * d[3] * (x[13] * x[13]) - e[3] + t9;
			con[4] = -2 * t5 - 3 * d[4] * (x[14] * x[14]) - e[4] + t10;

			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;

		}
	}
}

