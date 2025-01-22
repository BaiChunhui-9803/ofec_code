#include "f8_hybrid3.h"
#include "../../global/classical/scaffer_F6.h"
#include "../../global/classical/hg_bat.h"
#include "../../global/classical/rosenbrock.h"
#include "../../global/classical/modified_schwefel.h"
#include "../../global/classical/elliptic.h"

namespace ofec {
	namespace cec2015 {
		void F8_hybrid3::setFunction() {
			m_num_function = 5;
			m_function.resize(m_num_function);
			m_param_fun.resize(m_num_function);
			m_start.resize(m_num_function);
			m_dim.resize(m_num_function);
			m_hybrid_translation.resize(m_num_function);
			m_temp_sol.clear();

			size_t i, tmp;
			Real f_p[5] = { 0.1, 0.2, 0.2, 0.2, 0.3 };
			BasicFunctions f(5);
			f[0] = &createFunction<ScafferF6>;
			f[1] = &createFunction<HGBat>;
			f[2] = &createFunction<Rosenbrock>;
			f[3] = &createFunction<ModifiedSchwefel>;
			f[4] = &createFunction<Elliptic>;

			tmp = 0;
			for (i = 0; i < m_num_function - 1; ++i) {
				m_dim[i] = (size_t)ceil(f_p[i] * m_number_variables);
				tmp += m_dim[i];
			}
			m_dim[m_num_function - 1] = m_number_variables - tmp;
			m_start[0] = 0;
			for (i = 1; i < m_num_function; ++i)
				m_start[i] = m_start[i - 1] + m_dim[i - 1];
			for (i = 0; i < m_num_function; ++i)
				m_temp_sol.push_back(Solution<>(1, 0, m_dim[i]));
			for (size_t i = 0; i < m_num_function; ++i) {
				auto param = *m_param;
				param["problem name"] = m_name + "_part" + std::to_string(i + 1);
				param["number of variables"] = (int)m_dim[i];
				m_param_fun[i] = std::make_shared<const ParameterMap>(param);
				m_function[i].reset(dynamic_cast<Function*>(f[i]()));
				m_function[i]->setIdRnd(m_random.get());
				m_function[i]->setIdParam(m_param_fun[i]);
				m_function[i]->initialize();
				m_function[i]->setBias(0);
			}
		}

		void F8_hybrid3::evaluateObjective(Real *x, std::vector<Real> &obj) {
			Hybrid::evaluateObjective(x, obj);
			obj[0] += 800;
		}
	}
}