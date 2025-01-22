#include "f13_composition5.h"
#include "f8_hybrid3.h"
#include "../../global/classical/rastrigin.h"
#include "f6_hybrid1.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/scaffer_F6.h"

namespace ofec {
	namespace cec2015 {
		F13_global_composition5::F13_global_composition5(const ParameterMap &v) :
			F13_global_composition5((v.at("problem name")), (v.at("number of variables")), 1) {
			
		}
		F13_global_composition5::F13_global_composition5(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			composition_2015(name, size_var, size_obj) {
			
		}

		F13_global_composition5::~F13_global_composition5() {
			if (m_num_function > 0) {
				for (auto &i : m_hybrid)
					if (i) delete i;
			}
		}
		void F13_global_composition5::setFunction() {
			BasicFunctions f(5);
			f[0] = &createFunction<F8_hybrid3>;
			f[1] = &createFunction<rastrigin>;
			f[2] = &createFunction<F6_hybrid1>;
			f[3] = &createFunction<schwefel>;
			f[4] = &createFunction<scaffer_F6>;

			size_t num1 = 0;
			size_t num2 = 0;
			for (size_t i = 0; i < m_num_function + m_num_hybrid; ++i) {
				if(i==0||i==2) m_hybrid[num1++] = dynamic_cast<Hybrid*>(f[i]("", m_number_variables, m_number_objectives));
				else m_function[num2++] = dynamic_cast<function*>(f[i]("", m_number_variables, m_number_objectives));
				m_function[i]->initialize();
			}

			for (auto &i : m_function)
				i->setConditionNumber(1.);

			m_converge_severity[0] = 10;
			m_converge_severity[1] = 10;
			m_converge_severity[2] = 10;
			m_converge_severity[3] = 20;
			m_converge_severity[4] = 20;

			m_height[0] = 1;
			m_height[1] = 10;
			m_height[2] = 1;
			m_height[3] = 25;
			m_height[4] = 10;

			m_f_bias[0] = 0;
			m_f_bias[1] = 100;
			m_f_bias[2] = 200;
			m_f_bias[3] = 300;
			m_f_bias[4] = 400;

		}

		void F13_global_composition5::initialize() {
			m_variable_monitor = true;
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_function = 3;
			m_num_hybrid = 2;
			m_function.resize(m_num_function);
			m_hybrid.resize(m_num_hybrid);
			m_height.resize(m_num_function + m_num_hybrid);
			m_converge_severity.resize(m_num_function + m_num_hybrid);
			m_f_bias.resize(m_num_function + m_num_hybrid);
			setFunction();

			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			load_hybrid_translation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");

			for (auto &i : m_function) {
				i->get_optima().clear();
				i->setGlobalOpt(i->translation().data());
			}
			size_t num = 0;
			std::vector<Real> hy_bias = { 800,600 };
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			for (auto &i : m_hybrid) {
				for (size_t j = 0; j < m_number_variables; ++j)
					s.variable()[j] = i->get_optima().variable(0)[j] + i->hybrid_translation()[j];
				i->get_optima().clear();
				i->get_optima().append(s.variable());
				i->evaluate_(s, caller::Problem, false, false);
				s.objective()[0] -= hy_bias[num++];
				i->get_optima().append(s.objective());
			}
			// Set optimal solution

			s.variable() = m_hybrid[0]->get_optima().variable(0);
			m_optima->append(s.variable());

            s.evaluate(false, caller::Problem);
            m_optima->append(s.objective());

			add_tag(problem_tag::MMOP);
			m_initialized = true;
		}
		int F13_global_composition5::evaluateObjective(Real *x, std::vector<Real> &obj) {
			std::vector<Real> x_(m_number_variables);
			std::copy(x, x + m_number_variables, x_.begin());
			std::vector<Real> weight(m_num_function + m_num_hybrid, 0);

			set_weight(weight, x_);
			std::vector<Real> fit(m_num_function + m_num_hybrid);
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			size_t num1 = 0;
			size_t num2 = 0;
			std::vector<Real> hy_bias = { 800,600 };
			for (size_t i = 0; i < m_num_function + m_num_hybrid; ++i) { // calculate objective value for each function
				s.variable().vect() = x_;

				if (i == 0 || i == 2) {
					for (size_t j = 0; j < m_number_variables; ++j)
						s.variable()[j] -= m_hybrid[num1]->hybrid_translation()[j];
					m_hybrid[num1]->evaluate_(s, caller::Problem, false, false);
					fit[i] = s.objective()[0] - hy_bias[num1++];
				}
				else { 
					m_function[num2++]->evaluate_(s, caller::Problem, false, false); 
					fit[i] = s.objective()[0];
				}

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
			obj[0] += 1300;
			return kNormalEval;
		}
		bool F13_global_composition5::load_hybrid_translation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, m_name + "_HyShift_");
			s.insert(0, path);    // data path
			s.insert(0, g_working_directory);

			for (auto &i : m_hybrid)
				i->hybrid_translation().resize(m_number_variables);
			std::ifstream in(s);
			if (in.fail()) {
				setTranslation();
				std::ofstream out(s);
				for (size_t i = 0; i < m_num_hybrid; ++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						out << m_hybrid[i]->hybrid_translation()[j] << " ";

				out.close();
			}
			else {
				for (size_t i = 0; i < m_num_hybrid; ++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						in >> m_hybrid[i]->hybrid_translation()[j];
			}
			in.close();
			return true;
		}

		void F13_global_composition5::set_hybrid_translation() {
			for (int i = 0; i < m_num_hybrid; i++)
				for (int j = 0; j < m_number_variables; j++)
					m_hybrid[i]->hybrid_translation()[j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*(1 - global::ms_global->m_uniform[caller::Problem]->next());
		}
		void F13_global_composition5::set_weight(std::vector<Real>& weight, const std::vector<Real>&x) {
			size_t num1 = 0;
			size_t num2 = 0;
			for (size_t i = 0; i < m_num_function + m_num_hybrid; ++i) { // calculate weight for each function
				weight[i] = 0;
				for (size_t j = 0; j < m_number_variables; ++j) {
					//weight[i] += pow(x[j] - m_function[i]->translation()[j], 2);
					if(i==0||i==2) weight[i] += pow(x[j] - m_hybrid[num1]->get_optima().variable(0)[j], 2);
					else weight[i] += pow(x[j] - m_function[num2]->get_optima().variable(0)[j], 2);
				}
				if (i == 0 || i == 2) ++num1;
				else ++num2;
				if (fabs(weight[i])>1e-6) weight[i] = exp(-weight[i] / (2 * m_number_variables*m_converge_severity[i] * m_converge_severity[i])) / (Real)(pow(weight[i], 0.5));
				else {
					for (auto &m : weight) {
						m = 0;
					}
					weight[i] = 1;
					break;
				}
			}
		}
		Hybrid* F13_global_composition5::get_hybrid(size_t num) {
			return m_hybrid[num];
		}
	}
}