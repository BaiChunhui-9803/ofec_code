#include "F11_hybrid2.h"
#include "../../global/classical/griewank.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/rosenbrock.h"
#include "../../global/classical/scaffer_F6.h"
namespace ofec {
	namespace cec2015 {
		F11_hybrid2::F11_hybrid2(const ParameterMap &v) :
			F11_hybrid2((v.at("problem name")), (v.at("number of variables")), 1) {
			
		}
		F11_hybrid2::F11_hybrid2(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			hybrid(name, size_var, size_obj) {
			
		}
		void F11_hybrid2::setFunction() {
			size_t i, tmp;
			Real f_p[4] = { 0.2, 0.2, 0.3, 0.3 };
			BasicFunctions f(4);
			f[0] = &createFunction<griewank>;
			f[1] = &createFunction<weierstrass>;
			f[2] = &createFunction<rosenbrock>;
			f[3] = &createFunction<scaffer_F6>;

			tmp = 0;
			for (i = 0; i<m_num_function - 1; ++i)
			{
				m_dim[i] = (size_t)ceil(f_p[i] * m_number_variables);
				tmp += m_dim[i];
			}
			m_dim[m_num_function - 1] = m_number_variables - tmp;
			m_start[0] = 0;
			for (i = 1; i<m_num_function; ++i)
			{
				m_start[i] = m_start[i - 1] + m_dim[i - 1];
			}
			for (size_t i = 0; i < m_num_function; ++i) {
				m_function[i] = dynamic_cast<function*>(f[i]("", m_dim[i], m_number_objectives));
				m_function[i]->initialize();
				m_function[i]->setBias(0);
			}
		}

		void F11_hybrid2::initialize() {
			m_variable_monitor = true;
			m_num_function = 4;
			m_function.resize(m_num_function);
			m_start.resize(m_num_function);
			m_dim.resize(m_num_function);
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			//setBias(1100);
			setFunction();
			size_t count = 0;
			for (auto &i : m_random_perm)
				i = count++;
			global::ms_global->m_uniform[caller::Problem]->shuffle(m_random_perm.begin(), m_random_perm.end());
			// Set optimal solution
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			size_t num = 0;
			for (size_t i = 0; i < m_num_function;++i) {
				for (size_t j = 0; j < m_dim[i];++j)
                    s.variable()[m_random_perm[num++]] = m_function[i]->get_optima().variable(0)[j];
			}
			
			m_optima->append(s.variable());

			s.evaluate(false, caller::Problem);
			m_optima->append(s.objective());
			// end set
			m_initialized = true;
		}
		int F11_hybrid2::evaluateObjective(Real *x, std::vector<Real> &obj) {

			hybrid::evaluateObjective(x, obj);
			obj[0] += 1100;
			return kNormalEval;
		}
	}
}