#include "F12_hybrid3.h"
#include "../../global/classical/katsuura.h"
#include "../../global/classical/happy_cat.h"
#include "../../global/classical/griewank_rosenbrock.h"
#include "../../global/classical/modified_schwefel.h"
#include "../../global/classical/ackley.h"
namespace ofec {
	namespace cec2015 {
		F12_hybrid3::F12_hybrid3(const ParameterMap &v) :
			F12_hybrid3((v.at("problem name")), (v.at("number of variables")), 1) {
			
		}
		F12_hybrid3::F12_hybrid3(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			hybrid(name, size_var, size_obj) {
			
		}
		void F12_hybrid3::setFunction() {
			size_t i, tmp;
			Real f_p[5] = { 0.1, 0.2, 0.2, 0.2, 0.3 };
			BasicFunctions f(5);
			f[0] = &createFunction<katsuura>;
			f[1] = &createFunction<happy_cat>;
			f[2] = &createFunction<griewank_rosenbrock>;
			f[3] = &createFunction<modified_schwefel>;
			f[4] = &createFunction<ackley>;

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
		void F12_hybrid3::initialize() {
			m_variable_monitor = true;
			m_num_function = 5;
			m_function.resize(m_num_function);
			m_start.resize(m_num_function);
			m_dim.resize(m_num_function);
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			//setBias(1200);
			setFunction();
			size_t count = 0;
			for (auto &i : m_random_perm)
				i = count++;
			global::ms_global->m_uniform[caller::Problem]->shuffle(m_random_perm.begin(), m_random_perm.end());
			// Set optimal solution
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
            for (size_t i = 0; i < m_number_variables; ++i) {
                s.variable()[i] = 0;
            }

			m_optima->append(s.variable());

            s.evaluate(false, caller::Problem);
            m_optima->append(s.objective());
			// end set
			m_initialized = true;
		}
		int F12_hybrid3::evaluateObjective(Real *x, std::vector<Real> &obj) {

			hybrid::evaluateObjective(x, obj);
			obj[0] += 1200;
			return kNormalEval;
		}
	}
}