#include "F13_composition2015_C5.h"
#include "../../global/classical/rosenbrock.h"
#include "../../global/classical/HGBat.h"
#include "../../global/classical/rastrigin.h"
#include "../../global/classical/ackley.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/katsuura.h"
#include "../../global/classical/scaffer_F6.h"
#include "../../global/classical/griewank_rosenbrock.h"
#include "../../global/classical/happy_cat.h"
#include "../../global/classical/modified_schwefel.h"

namespace ofec {
	namespace cec2015 {
		F13_composition2015_C5::F13_composition2015_C5(const ParameterMap &v) :
			F13_composition2015_C5((v.at("problem name")), (v.at("number of variables")), 1) {
			
		}
		F13_composition2015_C5::F13_composition2015_C5(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			composition_2015(name, size_var, size_obj) {
			
		}

		void F13_composition2015_C5::setFunction() {
			BasicFunctions f(5);
			f[0] = &createFunction<rosenbrock>;
			f[1] = &createFunction<HGBat>;
			f[2] = &createFunction<rastrigin>;
			f[3] = &createFunction<ackley>;
			f[4] = &createFunction<weierstrass>;
			f[0] = &createFunction<katsuura>;
			f[1] = &createFunction<scaffer_F6>;
			f[2] = &createFunction<griewank_rosenbrock>;
			f[3] = &createFunction<happy_cat>;
			f[4] = &createFunction<modified_schwefel>;

			for (size_t i = 0; i < m_num_function; ++i) {
				m_function[i] = dynamic_cast<function*>(f[i / 2]("", m_number_variables, m_number_objectives));
				m_function[i]->setBias(0);
			}

			for (auto &i : m_function)
				i->setConditionNumber(1.);

			m_converge_severity[0] = 10;
			m_converge_severity[1] = 20;
			m_converge_severity[2] = 30;
			m_converge_severity[3] = 40;
			m_converge_severity[4] = 50;
			m_converge_severity[5] = 60;
			m_converge_severity[6] = 70;
			m_converge_severity[7] = 80;
			m_converge_severity[8] = 90;
			m_converge_severity[9] = 100;

			m_height[0] = 0.1;
			m_height[1] = 10;
			m_height[2] = 10;
			m_height[3] = 0.1;
			m_height[4] = 2.5;
			m_height[5] = 1e-3;
			m_height[6] = 100;
			m_height[7] = 2.5;
			m_height[8] = 10;
			m_height[9] = 1;

			m_function[0]->setScale(100. / 2.048);
			m_function[1]->setScale(20.);
			m_function[2]->setScale(100. / 5.12);
			m_function[3]->setScale(1.);
			m_function[4]->setScale(100./0.5);
			m_function[5]->setScale(20.);
			m_function[6]->setScale(1.);
			m_function[7]->setScale(20.);
			m_function[8]->setScale(20.);
			m_function[9]->setScale(1. / 10.);

			for (auto &i : m_f_bias)
				i = 0;

			//setBias(1300);
		}

		void F13_composition2015_C5::initialize() {
			m_variable_monitor = true;
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_function = 10;
			m_function.resize(m_num_function);
			m_height.resize(m_num_function);
			m_converge_severity.resize(m_num_function);
			m_f_bias.resize(m_num_function);
			setFunction();
			load_optima("/instance/problem/continuous/multi_modal/cec2015/");
			loadTranslation("/instance/problem/continuous/multi_modal/cec2015/");
			loadRotation("/instance/problem/continuous/multi_modal/cec2015/");
			for (auto &i : m_function)
				i->set_rotation_flag(false);
			for (auto &i : m_function)
				i->set_scale_flag(false);
			// 10 or 20 or 30 Dim : 10 gopt 

			evaluate_optima();

			add_tag(problem_tag::MMOP);
			m_initialized = true;
		}
		int F13_composition2015_C5::evaluateObjective(Real *x, std::vector<Real> &obj) {
			std::vector<Real> x_(m_number_variables);
			std::copy(x, x + m_number_variables, x_.begin());
			std::vector<Real> weight(m_num_function, 0);

			set_weight(weight, x_);
			std::vector<Real> fit(m_num_function);
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			for (size_t i = 0; i < m_num_function; ++i) { // calculate objective value for each function
				s.variable().vect() = x_;
				for (size_t j = 0; j < m_number_variables; ++j)
					s.variable()[j] -= m_function[i]->translation()[j];
				rotate(i, s.variable().data());
				scale(i, s.variable().data());
				if (i == 0)
					for (size_t j = 0; j < m_number_variables; ++j)
						s.variable().data()[j] += 1;
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
			obj[0] += 1300;
			return kNormalEval;
		}

		bool F13_composition2015_C5::load_optima(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, m_name + "_Optima_");

			s.insert(0, path);// data path
			s.insert(0, g_working_directory);

			load_optima_(s);

			return true;
		}

		void F13_composition2015_C5::load_optima_(const std::string &path) {
			std::ifstream in;
			in.open(path);
			if (in.fail()) {
				set_optima();
				std::ofstream out(path);
				for (size_t i = 0; i < m_optima->number_variable(); ++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						out << optima()->solution(i).variable()[j] << " ";
				out.close();
			}
			else {
				while (!(in.eof())) {
					VariableVector<Real> temp_var(m_number_variables);
					for (size_t j = 0; j < m_number_variables; ++j)
						in >> temp_var[j];
					m_optima->append(std::move(temp_var));
				}
			}
			in.close();
		}

		void F13_composition2015_C5::set_optima() {
			//  to do ..
			throw myexcept("Waiting to do !");
		}

		bool F13_composition2015_C5::loadTranslation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << ".txt";
			s = ss.str();
			s.insert(0, m_name + "_Shift");
			s.insert(0, path);    // data path
			s.insert(0, g_working_directory);

			for (auto &i : m_function)
				i->translation().resize(m_number_variables);
			std::ifstream in(s);
			if (in.fail()) {
				throw myexcept("open translation file failed !");
				/*setTranslation();
				std::ofstream out(s);
				for (size_t i = 0; i < m_num_function; ++i) {
				for (size_t j = 0; j < m_number_variables; ++j)
				out << m_function[i]->translation()[j] << " ";
				}

				out.close();*/
			}
			else {
				for (size_t i = 0; i < m_num_function; ++i) {
					size_t count = 0;
					Real temp;
					while (count < 100) {
						++count;
						if (count <= m_number_variables)
							in >> m_function[i]->translation()[count - 1];
						else {
							in >> temp;
						}
					}
				}
			}
			in.close();

			return true;
		}
		void F13_composition2015_C5::setTranslation() {
			for (int i = 0; i < m_num_function; i++)
				for (int j = 0; j < m_number_variables; j++)
					m_function[i]->translation()[j] = (global::ms_global->m_uniform[caller::Problem]->next() - 0.5) * 2 * 80.;
		}

		void F13_composition2015_C5::evaluate_optima() {
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
            for (size_t i = 0; i < m_optima->number_variable(); ++i) {
                s.variable() = optima()->solution(i).variable();
                s.evaluate(false, caller::Problem);
                m_optima->append(s.objective());
            }
		}
		void F13_composition2015_C5::rotate(size_t num, Real *x) {
			Real *x_ = new Real[m_number_variables];
			std::copy(x, x + m_number_variables, x_);

			for (size_t i = 0; i<m_number_variables; ++i) {
				x[i] = 0;

				for (size_t j = 0; j < m_number_variables; ++j) {
					x[i] += m_function[num]->rotation()[i][j] * x_[j];
				}
			}

			delete[] x_;
			x_ = 0;
		}

		void F13_composition2015_C5::scale(size_t num, Real *x) {
			for (size_t i = 0; i<m_number_variables; ++i)
				x[i] /= m_function[num]->scale();
		}
		void F13_composition2015_C5::set_weight(std::vector<Real>& weight, const std::vector<Real>&x) {

			for (size_t i = 0; i < m_num_function; ++i) { // calculate weight for each function
				weight[i] = 0;
				for (size_t j = 0; j < m_number_variables; ++j) {
					//weight[i] += pow(x[j] - m_function[i]->translation()[j], 2);
					weight[i] += pow(x[j] - m_function[i]->get_optima().variable(0)[j], 2);
				}

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
	}
}