#include "Multi_DOP.h"
#include<algorithm>
#include<fstream>
#include "../../../../../utility/functional.h"
//#include<utility>

namespace ofec {
	Multi_DOP::Multi_DOP(const ParameterMap &v) :Multi_DOP(v.at("problem name"), v.at("number of variables"), v.at("numObj"), v.at("testItems"), v.at("numPeak")) {

	}

	Multi_DOP::Multi_DOP(const std::string &name, size_t size_var, size_t size_obj, size_t type, size_t size_peak) : problem(name, size_var, size_obj), DMOPs(name, size_var, size_obj) {

	}

	void Multi_DOP::initialize() {
		for (size_t i = 0; i < m_number_variables; ++i) {//all variables are in [-1,1]
			m_domain.setDomain(-1., 1., i);
		}

		set_change_fre(global::ms_arg.find("changeFre")->second);
		set_change_severity(global::ms_arg.find("changeSeverity")->second);
		//num_pri = std::floor(m_number_variables/(Real)(m_number_objectives-1));
		m_num_pri = 1;
		m_num_pub = m_number_variables - m_num_pri * (global::ms_arg.find("numObj")->second - 1);
		//m_num_pri = m_number_variables - m_num_pub;
		//m_num_PS = 2;
		//m_num_PS = global::ms_arg.find("numPS")->second;

		//初始化H，W，m_variables，R，需要将最高峰值放在第一个位置,H=[1,10],W=[5,45],R=[0.1,0.5]，m_variables=[-1,1]
		for (size_t i = 0; i < m_number_objectives; i++) {
			std::vector<Real> H_temp, W_temp;
			std::vector<std::vector<Real>> P_temp;
			for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
				H_temp.push_back(1 + 9 * global::ms_global->m_uniform[caller::Problem]->next());
				auto it = std::max_element(H_temp.begin(), H_temp.end());
				Real temp = *it;
				*it = *(H_temp.begin());
				*(H_temp.begin()) = temp;//将最高峰放在第一个元素位置
				W_temp.push_back(40 * global::ms_global->m_uniform[caller::Problem]->next() + 5);
				std::vector<Real> X_temp;
				for (size_t k = 0; k < m_number_variables; k++) {
					if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= m_number_objectives*m_num_pri)
						X_temp.push_back(2 * global::ms_global->m_uniform[caller::Problem]->next() - 1);
					else
						X_temp.push_back(0);
				}
				P_temp.push_back(X_temp);
			}
			m_H0.push_back(H_temp);
			m_W0.push_back(W_temp);
			m_x0.push_back(P_temp);
		}

		m_H = m_H0;
		m_W = m_W0;
		m_variables = m_x0;
		//generateAdLoadPF();
		m_initialized = true;
	}

	void Multi_DOP::set_change_type(size_t i) {
		dynamic_factor temp = dynamic_factor(i);
		dynamic_type = temp;
	}

	void Multi_DOP::update_problem(Real t) {
		std::string temp = global::ms_arg.find("testItems")->second;//temp的格式为：“1_2_13”
		std::vector<std::string> s;//以字符串方式存放索引
		std::string s1;
		for (auto it = temp.begin(); it != temp.end(); it++) {
			if (*it != '_')
				s1.push_back(*it);
			else {
				s.push_back(s1);
				s1.clear();
			}
			if (it == temp.end() - 1)
				s.push_back(s1);
		}

		for (auto &i : s) {//依次索引测试项目的元素，
			set_change_type(atoi(i.c_str()));
			switch (get_dynamic_type()) {
			case dynamic_factor::Change_pattern://只变问题变化的模式:state-dependency
			{
				t = 1. / global::ms_arg.find("changeSeverity")->second*(std::exp(1 / 5) - 1);
				break;
			}
			case dynamic_factor::Change_num_objective://只变子目标的个数,子目标数目变化，可以考虑同时变化PS的位置和半径来模拟PS的退化
			{
				size_t M = global::ms_arg.find("numObj")->second + std::floor(t / 2.);
				resize_objective(M);
				break;
			}
			case dynamic_factor::Change_num_pri://只变私有变量的个数
			{
				m_num_pri = 4 + 3 * std::sin(0.5*OFEC_PI*t);
				resize_variable(m_num_pri + m_num_pub);
				break;
			}
			case dynamic_factor::Change_num_pub://只变公有变量的个数
			{
				m_num_pub = 4 + 3 * std::sin(0.5*OFEC_PI*t);
				//num_pub = m_number_objectives + std::floor((m_number_objectives - 1) * std::sin(0.5*OFEC_PI*t));
				resize_variable(m_num_pri + m_num_pub);
				break;
			}
			case dynamic_factor::Change_location_peak://只变子目标峰的位置，同时要更新PS中心的位置
			{
				for (size_t i = 0; i < m_number_objectives; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						for (size_t k = 0; k < m_number_variables; k++) {
							if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= m_number_objectives*m_num_pri)
								m_variables[i][j][k] = m_x0[i][j][k] * std::cos(OFEC_PI*t);
						}
					}
				}
				break;
			}
			case dynamic_factor::Change_width_peak://只变子目标峰的宽度
			{
				for (size_t i = 0; i < m_number_objectives; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						m_W[i][j] = m_W0[i][j] * std::cos(OFEC_PI*t);
					}
				}
				break;
			}
			case dynamic_factor::Change_detect://只变问题的可检测性
			{
				size_t temp = std::ceil(global::ms_arg.find("numPeak")->second *global::ms_global->m_uniform[caller::Problem]->next());
				for (size_t i = 0; i < m_H.size(); i++) {
					m_H[i][temp] = m_H[i][temp] * (1 + 0.1*std::sin(0.5*OFEC_PI*t));//这里有可能导致第一个峰不再是最高峰
				}
				break;
			}
			case dynamic_factor::Change_predict://只变问题的可预测性，PS的位置加随机，更新
			{
				for (size_t i = 0; i < m_number_objectives; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						for (size_t k = 0; k < m_number_variables; k++) {
							if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= m_number_objectives*m_num_pri)
								if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= m_number_objectives*m_num_pri)
									m_variables[i][j][k] = 0.5*m_x0[i][j][k] * std::cos(OFEC_PI*t) + 0.5*(2 * global::ms_global->m_uniform[caller::Problem]->next() - 1);
						}
					}
				}
				break;
			}
			default:
			{
				break;
			}
			}
		}
	}

	int Multi_DOP::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();

		if (time_changed() && t != 0. && (!get_updated_state())) {//防止不计数评价重复更新问题和重复采样PF
			update_problem(t);
			/*m_optima.reset(new Optima<>());
			generateAdLoadPF();*/
			set_updated_state(true);
		}
		else if (m_evaluations % (size_t)global::ms_arg.find("changeFre")->second != 0) {
			set_updated_state(false);
		}

		size_t p = global::ms_arg.find("numPeak")->second;
		std::vector<int> num;//store the number of variables of each subobjective

		std::vector<Real> opt_pri;//store the optimum of all private dimension
		std::vector<std::vector<Real>> opt_pub;//store the optimum of all public dimension

		//直接将每个子目标的第一个峰的位置当做最高峰
		for (size_t i = 0; i < m_number_objectives; i++) {
			Real temp = 0;
			std::vector<Real> temp1;

			//计算解的目标值
			for (size_t j = 0; j < p; j++) {
				Real temp2 = 0;
				for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
					temp2 += std::pow(x[k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
					if (j == 0)//默认第一个位置为最高峰
						opt_pri.push_back(m_variables[i][j][k]);
				}
				for (size_t m = m_num_pri * m_number_objectives; m < m_number_variables; m++) {//public variables
					temp2 += std::pow(x[m] - m_variables[i][j][m], 2);
					if (j == 0)//默认第一个位置为最高峰
						temp1.push_back(m_variables[i][j][m]);
				}
				temp = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(temp2)), temp);
			}

			obj[i] = (m_H[i][0] -temp) / m_H[i][0];//归一化
			//obj[i] = (m_H[i][0] - boundary_value - std::fabs(temp - boundary_value)) / (m_H[i][0] - boundary_value);//归一化至[0,1]
			if (!temp1.empty())
				opt_pub.push_back(temp1);
		}
		return kNormalEval;
		//是否验证一下子目标最高峰位置张成的区域属于最优解集的可能性有多少？
	}
}