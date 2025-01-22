#include "f12_composition4.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/rastrigin.h"
#include "../../global/classical/elliptic.h"
#include "../../global/classical/scaffer_F6.h"
#include "../../global/classical/happy_cat.h"

namespace ofec {
	namespace cec2015 {
		void F12_composition4::setFunction() {
			BasicFunctions f(5);
			f[0] = &createFunction<schwefel>;
			f[1] = &createFunction<rastrigin>;
			f[2] = &createFunction<elliptic>;
			f[3] = &createFunction<scaffer_F6>;
			f[4] = &createFunction<happy_cat>;

			for (size_t i = 0; i < m_num_function; ++i) {
				m_function[i] = dynamic_cast<function*>(f[i]("", m_number_variables, m_number_objectives));
				m_function[i]->initialize();
				m_function[i]->setBias(0);
			}

			for (auto &i : m_function)
				i->setConditionNumber(1.);

			m_converge_severity[0] = 10;
			m_converge_severity[1] = 20;
			m_converge_severity[2] = 20;
			m_converge_severity[3] = 30;
			m_converge_severity[4] = 30;

			m_height[0] = 0.25;
			m_height[1] = 1;
			m_height[2] = 1e-7;
			m_height[3] = 10;
			m_height[4] = 10;

			m_f_bias[0] = 0;
			m_f_bias[1] = 100;
			m_f_bias[2] = 100;
			m_f_bias[3] = 200;
			m_f_bias[4] = 200;

		}

		void F12_composition4::initialize() {
			m_num_function = 5;
			m_function.resize(m_num_function);
			m_height.resize(m_num_function);
			m_converge_severity.resize(m_num_function);
			m_f_bias.resize(m_num_function);
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_variable_monitor = true;
			setFunction();

			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
			for (auto &i : m_function) {
				i->get_optima().clear();
				i->setGlobalOpt(i->translation().data());
			}
			// Set optimal solution
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
            s.variable() = m_function[0]->get_optima().variable(0);

            m_optima->append(s.variable());
            m_optima->setVariableGiven(true);

            s.evaluate(false, caller::Problem);
            m_optima->append(s.objective());

			add_tag(problem_tag::MMOP);
			m_initialized = true;
		}
		int F12_composition4::evaluateObjective(Real *x, std::vector<Real> &obj) {
			std::vector<Real> x_(m_number_variables);
			std::copy(x, x + m_number_variables, x_.begin());
			std::vector<Real> weight(m_num_function, 0);

			set_weight(weight, x_);
			std::vector<Real> fit(m_num_function);
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			for (size_t i = 0; i < m_num_function; ++i) { // calculate objective value for each function
				s.variable().vect() = x_;

				m_function[i]->evaluate_(s, caller::Problem, false, false);
				fit[i] = s.objective()[0];

			}
			Real sumw = 0;
			for (size_t i = 0; i < m_num_function; ++i)
				sumw += weight[i];
			for (size_t i = 0; i < m_num_function; ++i)
				weight[i] /= sumw;

			Real temp = 0;
			for (size_t i = 0; i < m_num_function; ++i) {
				temp += weight[i] * (m_height[i] * fit[i] + m_f_bias[i]);
			}

			obj[0] = temp;
			obj[0] += 1200;
			return kNormalEval;
		}

	}
}