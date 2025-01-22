#include "SDMOP.h"
#include<algorithm>
#include<fstream>
#include "../../../../../utility/functional.h"
//#include<utility>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	SDMOP::SDMOP(const ParameterMap &v) :SDMOP(v.at("problem name"), v.at("number of variables"), v.at("numObj"), v.at("testItems"), v.at("numPeak")) {

	}

	SDMOP::SDMOP(const std::string &name, size_t size_var, size_t size_obj, size_t type, size_t size_peak) : problem(name, size_var, size_obj), DMOPs(name, size_var, size_obj) {

	}

	void SDMOP::initialize() {
		for (size_t i = 0; i < m_number_variables; ++i) {//all variables are in [-1,1]
			m_domain.setDomain(-1., 1., i);
		}

		set_change_fre(global::ms_arg.find("changeFre")->second);
		set_change_severity(global::ms_arg.find("changeSeverity")->second);
		//num_pri = std::floor(m_number_variables/(Real)(m_number_objectives-1));
		m_num_pri = 1;
		m_num_pub = m_number_variables-m_num_pri*(global::ms_arg.find("numObj")->second-1);
		//m_num_pri = m_number_variables - m_num_pub;
		//m_num_PS = 2;
		m_num_PS = global::ms_arg.find("numPS")->second;

		//��ʼ��H��W��m_variables��R����Ҫ����߷�ֵ���ڵ�һ��λ��,H=[1,10],W=[5,45],R=[0.1,0.5]��m_variables=[-1,1]
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			std::vector<Real> H_temp, W_temp;
			std::vector<std::vector<Real>> P_temp;
			for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
				H_temp.push_back(1 + 9 * global::ms_global->m_uniform[caller::Problem]->next());
				auto it = std::max_element(H_temp.begin(), H_temp.end());
				Real temp = *it;
				*it = *(H_temp.begin());
				*(H_temp.begin()) = temp;//����߷���ڵ�һ��Ԫ��λ��
				W_temp.push_back(40 * global::ms_global->m_uniform[caller::Problem]->next() + 5);
				std::vector<Real> X_temp;
				for (size_t k = 0; k < m_number_variables; k++) {
					if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= (m_number_objectives - 1)*m_num_pri)
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

		//��ʼ��PS�ĸ���������λ�ü��뾶
		for (size_t i = 0; i < m_num_PS; i++) {
			if (i < m_number_objectives - 1) {
				std::vector<Real> temp1;
				for (size_t j = 0; j < m_num_pub; j++) {
					//����Ŀ������ֵ��λ�ø�ֵ������
					temp1.push_back(m_x0[i][0][m_num_pri*(m_number_objectives - 1) + j]);
				}
				m_center_of_PS0.push_back(temp1);
			}
			else {
				std::vector<Real> temp2;
				for (size_t k = 0; k < m_num_pub; k++) {
					//�������ĵ�λ�������ֵ
					temp2.push_back(2 * global::ms_global->m_uniform[caller::Problem]->next() - 1);
				}
				m_center_of_PS0.push_back(temp2);
			}
			m_R0.push_back(global::ms_arg.find("radiusPS")->second);
			//m_R.push_back(0.9 / (1 + std::exp(-1 * (m_number_objectives - 3))));//PS�뾶��Ŀ�������

		}
		m_center_of_PS = m_center_of_PS0;
		m_R = m_R0;

		generateAdLoadPF();
		m_initialized = true;
	}

	void SDMOP::generateAdLoadPF() {
		/**********************************************************************
		��Ŀ��ÿһά�Ĳ�������num��Ŀ����M�Ĺ�ϵͨ��������϶��ɣ�(�ܹ���Ҫ��M-1ά)
		1����M=2ʱ��num=2000;   total_num=2000;
		2����2<M<8ʱ��num=ceil(-3.333*M.^3+58.93*M.^2-348.5*M+704.9);  total_num=num^(M-1);
		3����M>=8ʱ��num=10;     total_num=num^(M-1);
		***********************************************************************/
		//std::vector<std::vector<Real>> sample_PF = get_1d_sample_PF();
		//Real t = get_time();
		//std::vector<std::vector<Real>> sample_PF = get_1d_sample_PF(m_number_objectives);
		//matrix o_point(sample_PF[0].size(), m_number_objectives);
		Real m_offset = std::min(1 - m_T / 2. * sign(2.* std::sqrt(m_number_objectives-1) / m_T - 2 * std::ceil(std::sqrt(m_number_objectives - 1) / m_T) + 1)*std::pow(std::fabs(2 * std::sqrt(m_number_objectives - 1) / m_T - 2. * std::ceil(std::sqrt(m_number_objectives - 1) / m_T) + 1), m_alpha) - m_T * (std::ceil(std::sqrt(m_number_objectives - 1) / m_T) - 1. / 2), 0.);
		size_t num;
		if (m_number_objectives == 2) {
			//std::vector<std::vector<Real>> sam_PF = get_1d_sample_PF(2000);
			get_1d_sample_PF(2000);
			matrix o_point(m_sample_PF[0].size(), m_number_objectives);
			//for (int i = 0; i < m_sample_PF[0].size(); i++) {
			//	o_point[i][0] = m_sample_PF[0][i];
			//	o_point[i][1] = 1 - m_T / 2. * sign(2.* o_point[i][0] / m_T - 2 * std::ceil(o_point[i][0] / m_T) + 1)*std::pow(std::fabs(2 * o_point[i][0] / m_T - 2. * std::ceil(o_point[i][0] / m_T) + 1), m_alpha) - m_T * (std::ceil(o_point[i][0] / m_T) - 1. / 2) - m_offset;
			//	m_optima->append(o_point[i].vect());
			//	m_optima->set_flag_objective(true);//given optima
			//}
			std::ofstream out("G:/console/result/PF.txt");
			if (out) {
			    out << "f1" << " " << "f2" << std::endl;
			    for (int i = 0; i < m_sample_PF[0].size(); i++) {
			        o_point[i][0] = m_sample_PF[0][i];
			        o_point[i][1] = 1 - m_T / 2. * sign(2.* o_point[i][0] / m_T - 2 * std::ceil(o_point[i][0] / m_T) + 1)*std::pow(std::fabs(2 * o_point[i][0] / m_T - 2. * std::ceil(o_point[i][0] / m_T) + 1), m_alpha) - m_T * (std::ceil(o_point[i][0] / m_T) - 1. / 2) - m_offset;
			        m_optima->append(o_point[i].vect());
					m_optima->set_flag_objective(true);//given optima
			        out << o_point[i][0] << " " << o_point[i][1] << std::endl;
			     }
			     out.close();
			}
		}
		else {
			if (m_number_objectives >= 8) {
				num = 10;
			}
			else {
				num = std::ceil(-3.333*std::pow(m_number_objectives, 3) + 58.93*std::pow(m_number_objectives, 2) - 348.5*m_number_objectives + 704.9);
			}
			//std::vector<std::vector<Real>> sample_PF = get_1d_sample_PF(num);
			get_1d_sample_PF(num);
			size_t sum = 1;//�ܵĲ�������
			for (auto &i : m_sample_PF) {
				sum *= i.size();
			}
			matrix o_point(sum, m_number_objectives);
			//for (size_t j = 0; j < sum; j++) {
			//	Real temp = 0;
			//	for (size_t k = 0; k < m_number_objectives - 1; k++) {
			//		size_t product = 1;
			//		for (size_t m = k + 1; m < m_number_objectives - 1; m++) {
			//			product *=m_sample_PF[m].size();
			//		}
			//		size_t ind;
			//		if (k < m_number_objectives - 2) {
			//			ind = j / product;
			//		}
			//		else {
			//			ind = j % m_sample_PF[k].size();
			//		}
			//		o_point[j][k] = m_sample_PF[k][ind];
			//		temp += std::pow(o_point[j][k], 2);
			//	}
			//	temp = std::sqrt(temp);//������ת�õ��ĵ�Чֵ
			//	o_point[j][m_number_objectives - 1] = 1 - m_T / 2. * sign(2.* temp / m_T - 2 * std::ceil(temp / m_T) + 1)*std::pow(std::fabs(2 * temp / m_T - 2. * std::ceil(temp / m_T) + 1), m_alpha) - m_T * (std::ceil(temp / m_T) - 1. / 2) - m_offset;;
			//	m_optima->append(o_point[j].vect());
			//	m_optima->set_flag_objective(true);//given optima
			//}
			std::ofstream out("G:/console/result/PF.txt");
			if (out) {
			    out << "f1" << " " << "f2" << " "<<"f3"<<std::endl;
			    for (size_t j = 0; j < sum; j++) {
			        Real temp = 0;
			        for (size_t k = 0; k < m_number_objectives - 1; k++) {
			            size_t product = 1;
			            for (size_t m = k + 1; m < m_number_objectives - 1; m++) {
			                product *=m_sample_PF[m].size();
			            }
			            size_t ind;
			            if (k < m_number_objectives - 2) {
			               ind = j / product;
			            }
			            else {
							ind = j % m_sample_PF[k].size();
			            }
			            o_point[j][k] = m_sample_PF[k][ind];
			            temp += std::pow(o_point[j][k], 2);
			        }
			        temp = std::sqrt(temp);//������ת�õ��ĵ�Чֵ
			        o_point[j][m_number_objectives - 1] = 1 - m_T / 2. * sign(2.* temp / m_T - 2 * std::ceil(temp / m_T) + 1)*std::pow(std::fabs(2 * temp / m_T - 2. * std::ceil(temp / m_T) + 1), m_alpha) - m_T * (std::ceil(temp / m_T) - 1. / 2) - m_offset;;
                    m_optima->append(o_point[j].vect());
					m_optima->set_flag_objective(true);//given optima
			        out << o_point[j][0] << " " << o_point[j][1] <<  " " << o_point[j][2] <<std::endl;
			    }
			    out.close();
			}
		}
	}

	void SDMOP::get_1d_sample_PF(size_t n) {//n is the number of sample points in each objective dimension
		size_t p = global::ms_arg.find("numPeak")->second;
		m_PS_obj_range.clear();
		m_boundary_point.clear();
		//*****�������Ϸ��۵��ٽ�ֵ
		std::vector<Real> m_obj_boundary_value;//�����Ŀ�겿�������߽�˵��ƽ��ֵ
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			Real boundary_value = 0;//��Ų��������߽�˵��ƽ��ֵ
		    //��ȡ�߽������
			std::vector<std::vector<Real>> boundary_point;
			for (size_t q = 0; q < 2 * m_number_variables; q++) {
				std::vector<Real> point;
				for (size_t n = 0; n < m_number_variables; n++) {
					Real tmp3;
					if (global::ms_global->m_uniform[caller::Problem]->next() > 0.5)
						tmp3 = 1.;
					else
						tmp3 = -1.;
					point.push_back(tmp3);
				}
				boundary_point.push_back(point);
			}
			if (0) {
				//����߽綥�����Сֵ
				boundary_value = m_H[i][0];
				for (size_t q = 0; q < boundary_point.size(); q++) {
					Real tmp1 = 0;
					for (size_t j = 0; j < p; j++) {
						Real tmp2 = 0;
						for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
							tmp2 += std::pow(boundary_point[q][k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
						}
						for (size_t m = m_num_pri * (m_number_objectives - 1); m < m_number_variables; m++) {//public variables
							tmp2 += std::pow(boundary_point[q][m] - m_variables[i][j][m], 2);
						}
						tmp1 = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(tmp2)), tmp1);
					}
					if (boundary_value > tmp1)
						boundary_value = tmp1;
				}
				m_obj_boundary_value.push_back(boundary_value);
			}
			else {
				//����߽������ƽ��ֵ
				for (size_t q = 0; q < boundary_point.size(); q++) {
					Real tmp1 = 0;
					for (size_t j = 0; j < p; j++) {
						Real tmp2 = 0;
						for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
							tmp2 += std::pow(boundary_point[q][k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
						}
						for (size_t m = m_num_pri * (m_number_objectives - 1); m < m_number_variables; m++) {//public variables
							tmp2 += std::pow(boundary_point[q][m] - m_variables[i][j][m], 2);
						}
						tmp1 = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(tmp2)), tmp1);
					}
					boundary_value += tmp1;
				}
				boundary_value = boundary_value / boundary_point.size();
				m_obj_boundary_value.push_back(boundary_value);
			}
		}

        //�ȼ���ÿ����Ŀ����ÿ��PS�������С���ֵ
		Obj_range obj_boundary_points;
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			std::vector<std::pair<Real, Real>> limit_value;
			for (size_t j = 0; j < m_num_PS; j++) {
				Real max_max = 0;
				Real max_min = 0;
				Real m_max, m_min;
				for (size_t k = 0; k < p; k++) {
					Real temp = 0;
					for (size_t m = 0; m < m_num_pub; m++) {
						temp += std::pow(m_variables[i][k][m_num_pri*(m_number_objectives - 1) + m] - m_center_of_PS[j][m], 2);
					}
					m_min = m_H[i][k] / (1 + m_W[i][k] * (std::sqrt(temp) + m_R[j]));//��j��PS�����ڵ�k�����µ���Сֵ,��Ĭ�ϵ�һ�������
					m_max = m_H[i][k] / (1 + m_W[i][k] * std::max(std::sqrt(temp) - m_R[j], (float)0.));//��j��PS�����ڵ�k�����µ����ֵ
					if (m_max > max_max)
						max_max = m_max;//ѡ��ÿ�����µ�������޺�������Ϊ���PS��ȡֵ��Χ
					if (m_min > max_min)
						max_min = m_min;
				}
				auto temp1 = std::make_pair(max_min, max_max);
				limit_value.push_back(temp1);//δ��һ��
			}
			obj_boundary_points.push_back(limit_value);
			limit_value.clear();

			//�洢��ת��ķ�Χ��
			for (size_t j = 0; j < m_num_PS; j++) {
				Real temp;
				if (obj_boundary_points[i][j].first < m_obj_boundary_value[i])  //�ж���Ŀ��ȡֵ�����Ƿ�С�ڱ߽緭�۵��ֵ
					if (obj_boundary_points[i][j].second < m_obj_boundary_value[i]) {
						temp = obj_boundary_points[i][j].first;
						obj_boundary_points[i][j].first = m_obj_boundary_value[i] - obj_boundary_points[i][j].second;
						obj_boundary_points[i][j].second= m_obj_boundary_value[i] - temp;
					}
					else if (m_obj_boundary_value[i] - obj_boundary_points[i][j].first > obj_boundary_points[i][j].second- m_obj_boundary_value[i]) {
						obj_boundary_points[i][j].second = m_obj_boundary_value[i] - obj_boundary_points[i][j].first;
						obj_boundary_points[i][j].first = 0;
					}
					else {
						obj_boundary_points[i][j].first = 0;
						obj_boundary_points[i][j].second = obj_boundary_points[i][j].second - m_obj_boundary_value[i];
					}
				else {
					obj_boundary_points[i][j].first = obj_boundary_points[i][j].first - m_obj_boundary_value[i];
					obj_boundary_points[i][j].second = obj_boundary_points[i][j].second - m_obj_boundary_value[i];
				}

				//����Ŀ��PS�����ȡֵ��Χ�Ķ˵����
				auto temp1=std::make_pair(1-obj_boundary_points[i][j].second/(m_H[i][0]-m_obj_boundary_value[i]), 1 - obj_boundary_points[i][j].first / (m_H[i][0] - m_obj_boundary_value[i]));
				limit_value.push_back(temp1);
			}
			m_PS_obj_range.push_back(limit_value);//��һ��������

			//��������PS��ȡֵ��Χ���е�������
			for (size_t n = 0; n < m_num_PS - 1; n++) {
				std::pair<Real, Real> temp;
				for (size_t p = 0; p < m_num_PS - 1 - n; p++) {
					if (obj_boundary_points[i][p].first>obj_boundary_points[i][p+1].first) {
						temp = obj_boundary_points[i][p];
						obj_boundary_points[i][p] = obj_boundary_points[i][p+1];
						obj_boundary_points[i][p + 1] = temp;
					}
				}
			}
		}

		//��������֮���ÿ����Ŀ���ȡֵ����,���߽�㱣��
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			Real min_min = obj_boundary_points[i][0].first;//�������PS��������С������ֵ
			Real temp1 = obj_boundary_points[i][0].second;//�������ֵ��С��PS����Ӧ������ֵ
			std::vector<Real> boundary_point;//����Ƭ�ε������յ�ֵ��Ԫ�ظ���Ϊż��
			boundary_point.push_back(1-min_min / (m_H[i][0] - m_obj_boundary_value[i]));
			for (size_t j = 0; j < m_num_PS; j++) {
				if (obj_boundary_points[i][j].second > temp1) {
					if (obj_boundary_points[i][j].first > temp1) {
						boundary_point.push_back(1-temp1  / (m_H[i][0] - m_obj_boundary_value[i]));//�����һ���ı߽��.��С���Ż�
						boundary_point.push_back(1-obj_boundary_points[i][j].first / (m_H[i][0] - m_obj_boundary_value[i]));
						min_min = obj_boundary_points[i][j].first;
						temp1 = obj_boundary_points[i][j].second;
					}
					else {
						temp1 = obj_boundary_points[i][j].second;
					}
				}
				if (j == m_num_PS - 1)
					boundary_point.push_back(1-temp1 / (m_H[i][0] - m_obj_boundary_value[i]));
			}
			m_boundary_point.push_back(boundary_point);
		}

		//������Ŀ�������PS�ڵ�ȡֵ���
		m_sample_PF.clear();//�������M-1����Ŀ��Ĳ�����
		std::vector<Real> m_sum;//���ÿһ����Ŀ��ȡֵ�Ŀ��
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			Real sum = 0;
			for (size_t j = 0; j < m_boundary_point[i].size() - 1;) {
				sum += m_boundary_point[i][j] - m_boundary_point[i][j+1];
				j = j + 2;
			}
			m_sum.push_back(sum);
		}
		//������Ŀ��ֵ�Ŀ�ȴ�С����ÿ����Ŀ���ϵĲ��������
		std::vector<size_t> m_num_1d_obj;//���ÿһ����Ŀ���ϵĲ��������
		Real m_average=0;
		for (auto &i : m_sum) {
			m_average = m_average + i;
		}
		m_average = m_average / (m_number_objectives - 1);
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			size_t num = 2 + std::ceil(m_sum[i] / m_average * n);//ÿһά�����������˵�
			m_num_1d_obj.push_back(num);
		}
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			std::vector<Real> m_point;//���ĳһ��Ŀ���ϵĲ�����
			Real m_step = m_sum[i] / (m_num_1d_obj[i]-1);//����ᵼ�²���������΢�仯
			Real temp2;
			//��i����Ŀ�����
			for (int k = m_boundary_point[i].size() - 1; k > 0;) {
				temp2 = m_boundary_point[i][k];//��С���Ż�
				while (temp2 <= m_boundary_point[i][k - 1]) {
					m_point.push_back(temp2);
					temp2 += m_step;
				}
				m_point.push_back(m_boundary_point[i][k - 1]);//����ÿһ�εĶ˵�
				k = k - 2;
			}
			m_sample_PF.push_back(m_point);
		}
	}

	void SDMOP::set_change_type(size_t i) {
		dynamic_factor temp = dynamic_factor(i);
		dynamic_type = temp;
	}

	bool SDMOP::in_range(Real val, std::pair<Real, Real> range) {
		if (val <= range.second && val >= range.first)
			return true;
		else
			return false;
	}

	void SDMOP::update_problem(Real t) {
		std::string temp = global::ms_arg.find("testItems")->second;//temp�ĸ�ʽΪ����1_2_13��
		std::vector<std::string> s;//���ַ�����ʽ�������
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

		for (auto &i : s) {//��������������Ŀ��Ԫ�أ�
			set_change_type(atoi(i.c_str()));
			switch (get_dynamic_type()) {
			case dynamic_factor::Change_pattern://ֻ������仯��ģʽ:state-dependency
			{
				t = 1. / global::ms_arg.find("changeSeverity")->second*(std::exp(1 / 5) - 1);
				break;
			}
			case dynamic_factor::Change_PF_curvature://ֻ�仯PF������
			{
				m_alpha = std::pow(5, -1 * std::cos(0.5*OFEC_PI*t));
				break;
			}
			case dynamic_factor::Change_PF_seg://ֻ��PF�İ�͹����
			{
				m_T = std::pow(2, 3. / 2 * (std::sin(0.5*OFEC_PI*(t-1)) + 1) - 2);
				break;
			}
			case dynamic_factor::Change_num_objective://ֻ����Ŀ��ĸ���,��Ŀ����Ŀ�仯�����Կ���ͬʱ�仯PS��λ�úͰ뾶��ģ��PS���˻�
			{
				size_t M = global::ms_arg.find("numObj")->second + std::floor(t / 2.);
				resize_objective(M);
				break;
			}
			case dynamic_factor::Change_num_pri://ֻ��˽�б����ĸ���
			{
				m_num_pri = 4 + 3 * std::sin(0.5*OFEC_PI*t);
				resize_variable(m_num_pri + m_num_pub);
				break;
			}
			case dynamic_factor::Change_num_pub://ֻ�乫�б����ĸ���
			{
				m_num_pub = 4 + 3 * std::sin(0.5*OFEC_PI*t);
				//num_pub = m_number_objectives + std::floor((m_number_objectives - 1) * std::sin(0.5*OFEC_PI*t));
				resize_variable(m_num_pri + m_num_pub);
				break;
			}
			case dynamic_factor::Change_location_peak://ֻ����Ŀ����λ�ã�ͬʱҪ����PS���ĵ�λ��
			{
				for (size_t i = 0; i < m_number_objectives - 1; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						for (size_t k = 0; k < m_number_variables; k++) {
							if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= (m_number_objectives - 1)*m_num_pri)
								m_variables[i][j][k] = m_x0[i][j][k] * std::cos(OFEC_PI*t);
						}
					}
					//����PS���ĵ�λ��
					for (size_t k = 0; k < m_num_PS; k++) {
						for (size_t m = 0; m < m_num_pub; m++) {
							if (k < m_number_objectives - 1) {
								//����Ŀ������ֵ��λ�ø�ֵ������
								m_center_of_PS[k][m] = m_variables[i][0][m_num_pri*(m_number_objectives - 1) + m] + m_R[k] * std::sin(0.5*OFEC_PI*t);
							}
							else {
								//������������Ŷ�
								m_center_of_PS[k][m] = 0.5* m_center_of_PS0[k][m] + global::ms_global->m_uniform[caller::Problem]->next() - 0.5;
							}
						}
					}
					//R.push_back(0.5 / (1 + std::exp(-1 * (m_number_objectives - 3))));//PS�뾶��Ŀ�������
				}
				break;
			}
			case dynamic_factor::Change_width_peak://ֻ����Ŀ���Ŀ��
			{
				for (size_t i = 0; i < m_number_objectives - 1; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						m_W[i][j] = m_W0[i][j] * std::cos(OFEC_PI*t);
					}
				}
				break;
			}
			case dynamic_factor::Change_num_PS://ֻ��PS�ĸ�����ͬʱ��Ҫ����PS������λ�ú�PS�İ뾶
			{
				size_t num_PS_t = 1 + std::ceil((m_number_objectives - 1)*std::fabs(std::sin(0.1*OFEC_PI*t)));
				size_t temp = num_PS_t - m_num_PS;
				if (temp>0) {
					for (size_t i = 0; i < temp; i++) {
						//�������ĵ�λ�������ֵ
						std::vector<Real> temp1;
						for (size_t j = 0; j < m_num_pub; j++) {
							temp1.push_back(2 * global::ms_global->m_uniform[caller::Problem]->next() - 1);
						}
						m_center_of_PS.push_back(temp1);
					}
					m_R.push_back(0.9 / (1 + std::exp(-1 * (m_number_objectives - 3))));
				}
				else if (temp<0) {
					for (size_t j = 0; j < -1 * temp; j++) {
						m_center_of_PS.pop_back();
						m_R.pop_back();
					}
				}
				m_num_PS = num_PS_t;
				break;
			}
			case dynamic_factor::Change_center_PS://ֻ��PS������λ��
			{
				for (size_t i = 0; i < m_num_PS; i++) {
					for (size_t j = 0; j < m_num_pub; j++) {
						m_center_of_PS[i][j] = m_center_of_PS0[i][j] + m_R[i] * std::sin(OFEC_PI*t);//���ĵ�λ���ڰ뾶��Χ���Ŷ���ʹ����ֵ��λ��ʼ�ձ���
					}
				}
				break;
			}
			case dynamic_factor::Change_radius_PS://ֻ��PS�ĸ��ǰ뾶
			{
				for (size_t i = 0; i < m_R.size(); i++) {
					m_R[i] = 0.9 / (1 + std::exp(-1 * (m_number_objectives - 3))) + 0.9 / (1 + std::exp(-1 * (m_number_objectives - 3))) / 2 * std::sin(2 * OFEC_PI*t);
				}
				break;
			}
			case dynamic_factor::Change_detect://ֻ������Ŀɼ����
			{
				size_t temp = std::ceil(global::ms_arg.find("numPeak")->second *global::ms_global->m_uniform[caller::Problem]->next());
				for (size_t i = 0; i < m_H.size(); i++) {
					m_H[i][temp] = m_H[i][temp] * (1 + 0.1*std::sin(0.5*OFEC_PI*t));//�����п��ܵ��µ�һ���岻������߷�
				}
				break;
			}
			case dynamic_factor::Change_predict://ֻ������Ŀ�Ԥ���ԣ�PS��λ�ü����������
			{
				for (size_t i = 0; i < m_number_objectives; i++) {
					for (size_t j = 0; j < (size_t)global::ms_arg.find("numPeak")->second; j++) {
						for (size_t k = 0; k < m_number_variables; k++) {
							if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= (m_number_objectives - 1)*m_num_pri)
							if (i*m_num_pri <= k && k < (i + 1)*m_num_pri || k >= (m_number_objectives - 1)*m_num_pri)
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

	int SDMOP::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real t = get_time();

		if (time_changed() && t != 0. && (!get_updated_state())) {//��ֹ�����������ظ�����������ظ�����PF
			update_problem(t);
			m_optima.reset(new Optima<>());
			generateAdLoadPF();
			set_updated_state(true);
//#ifdef OFEC_DEMO
//			ofec_demo::g_buffer->appendProBuffer(this);
//#endif
		}
		else if (m_evaluations % (size_t)global::ms_arg.find("changeFre")->second != 0) {
			set_updated_state(false);
		}

		size_t p = global::ms_arg.find("numPeak")->second;
		std::vector<int> num;//store the number of variables of each subobjective

		std::vector<Real> opt_pri;//store the optimum of all private dimension
		std::vector<std::vector<Real>> opt_pub;//store the optimum of all public dimension

		//ֱ�ӽ�ÿ����Ŀ��ĵ�һ�����λ�õ�����߷�
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			Real temp = 0;
			std::vector<Real> temp1;

			//��ȡ�߽������
			std::vector<std::vector<Real>> boundary_point;
			for (size_t q = 0; q < 2 * m_number_variables; q++) {
				std::vector<Real> point;
				for (size_t n = 0; n < m_number_variables; n++) {
					Real tmp3;
					if (global::ms_global->m_uniform[caller::Problem]->next() > 0.5)
						tmp3 = 1.;
					else
						tmp3 = -1.;
					point.push_back(tmp3);
				}
				boundary_point.push_back(point);
			}
			Real boundary_value = 0;//��Ų��������߽�˵��ƽ��ֵ����Сֵ
			if (0) {
				//����߽綥�����Сֵ
				boundary_value = m_H[i][0];
				for (size_t q = 0; q < boundary_point.size(); q++) {
					Real tmp1 = 0;
					for (size_t j = 0; j < p; j++) {
						Real tmp2 = 0;
						for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
							tmp2 += std::pow(boundary_point[q][k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
						}
						for (size_t m = m_num_pri * (m_number_objectives - 1); m < m_number_variables; m++) {//public variables
							tmp2 += std::pow(boundary_point[q][m] - m_variables[i][j][m], 2);
						}
						tmp1 = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(tmp2)), tmp1);
					}
					if (boundary_value > tmp1)
						boundary_value = tmp1;
				}
			}
			else {
				//����߽������ƽ��ֵ
				for (size_t q = 0; q < boundary_point.size(); q++) {
					Real tmp1 = 0;
					for (size_t j = 0; j < p; j++) {
						Real tmp2 = 0;
						for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
							tmp2 += std::pow(boundary_point[q][k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
						}
						for (size_t m = m_num_pri * (m_number_objectives - 1); m < m_number_variables; m++) {//public variables
							tmp2 += std::pow(boundary_point[q][m] - m_variables[i][j][m], 2);
						}
						tmp1 = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(tmp2)), tmp1);
					}
					boundary_value += tmp1;
				}
				boundary_value = boundary_value / boundary_point.size();
			}

			//������Ŀ��ֵ
			for (size_t j = 0; j < p; j++) {
				Real temp2 = 0;
				for (size_t k = m_num_pri * i; k < m_num_pri*(i + 1); k++) {//private variables
					temp2 += std::pow(x[k] - m_variables[i][j][k], 2);//m_variables is the optimum in this dimension
					if (j == 0)//Ĭ�ϵ�һ��λ��Ϊ��߷�
						opt_pri.push_back(m_variables[i][j][k]);
				}
				for (size_t m = m_num_pri * (m_number_objectives - 1); m < m_number_variables; m++) {//public variables
					temp2 += std::pow(x[m] - m_variables[i][j][m], 2);
					if (j == 0)//Ĭ�ϵ�һ��λ��Ϊ��߷�
						temp1.push_back(m_variables[i][j][m]);
				}
				temp = std::max(m_H[i][j] / (1 + m_W[i][j] * std::sqrt(temp2)), temp);
			}
				
			//obj[i] = (m_H[i][0] -temp) / m_H[i][0];//��һ��
			obj[i] = 1-std::fabs(temp - boundary_value)/(m_H[i][0] - boundary_value);//��һ����[0,1]
			if (!temp1.empty())
				opt_pub.push_back(temp1);
		}

		Real F = 0;
		for (size_t k = 0; k < m_number_objectives - 1; k++) {
			F = F + std::pow(obj[k], 2);
		}
		F = std::sqrt(F);

		Real g = 1;
		for (size_t i = 0; i < opt_pri.size(); i++) {//private variables
			g += std::pow(x[i] - opt_pri[i], 2);
		}
		//g = std::sqrt(g);
		//Real temp3 = 1;
		std::vector<Real> temp_dist;//��ŵ㵽��PS����ľ���
		for (size_t j = 0; j < m_num_PS; j++) {//public variables
			Real temp = 0;
			for (size_t k = 0; k < m_num_pub; k++) {
				temp += std::pow(x[m_num_pri*(m_number_objectives - 1) + k] - m_center_of_PS[j][k], 2);
			}
			temp = std::max(std::sqrt(temp) - m_R[j],(float) 0.);
			temp_dist.push_back(temp);
			//temp3 *= temp;
		}
		Real temp3=temp_dist.back();
		for (auto &i : temp_dist) {
			if (i < temp3)//temp3Ϊ�����PS����С����
				temp3 = i;
		}
		g =  std::sqrt(g + temp3*temp3);
		/*if (g == 0.)
			g = g + 1;
		else
			g = g + 1 +2/(1+std::pow(2,-1*g));*/

		Real m_offset = std::min(1 - m_T / 2. * sign(2. * std::sqrt(m_number_objectives-1) / m_T - 2 * std::ceil(std::sqrt(m_number_objectives - 1) / m_T) + 1)*std::pow(std::fabs(2 * std::sqrt(m_number_objectives - 1) / m_T - 2. * std::ceil(std::sqrt(m_number_objectives - 1) / m_T) + 1), m_alpha) - m_T * (std::ceil(std::sqrt(m_number_objectives - 1) / m_T) - 1. / 2), 0.);
		Real m_F= m_T / 2. * sign(2.* F / m_T - 2 * std::ceil(F / m_T) + 1)*std::pow(std::fabs(2 * F / m_T - 2. * std::ceil(F / m_T) + 1), m_alpha) + m_T * (std::ceil(F / m_T) - 1. / 2) ;
		obj[m_number_objectives - 1] = g *(1 - m_F - m_offset);
		//obj[m_number_objectives - 1] = g * (1 - m_T / 2. * sign(2.* F / m_T - 2 * std::ceil(F / m_T) + 1)*std::pow(std::fabs(2 * F / m_T - 2. * std::ceil(F / m_T) + 1), m_alpha) - m_T * (std::ceil(F / m_T) - 1. / 2) - m_offset);

		//������PS�����ڵķ�֧���
		//���ж��Ƿ���ȡֵ������
		bool temp = 0;
		for (size_t i = 0; i < m_number_objectives - 1; i++) {
			for (size_t j = 0; j < m_num_PS; j++) {
				temp = in_range(obj[i], m_PS_obj_range[i][j]);
				if (temp)
					break;
			}
			if (!temp)
				break;
		}
		
		if (!temp) {
			//�Բ��������ڵ�Ŀ��ֵ��������
			Real temp_F = 0;
			for (size_t i = 0; i < m_number_objectives - 1; i++) {
				Real temp = 0;
				for (size_t j = 0; j < m_boundary_point[i].size(); j++) {
					if (obj[i] >= m_boundary_point[i][m_boundary_point[i].size() - j - 1])
						temp = m_boundary_point[i][m_boundary_point[i].size() - j - 1];
				}
				temp_F += std::pow(temp, 2);
			}
			temp_F = std::sqrt(temp_F);
			Real temp_obj= 1 - m_T / 2. * sign(2.* temp_F / m_T - 2 * std::ceil(temp_F / m_T) + 1)*std::pow(std::fabs(2 * temp_F / m_T - 2. * std::ceil(temp_F / m_T) + 1), m_alpha) - m_T * (std::ceil(temp_F / m_T) - 1. / 2) - m_offset;
			obj[m_number_objectives - 1] = 2 * temp_obj - obj[m_number_objectives - 1]/g;
			//obj[m_number_objectives - 1]=g*temp_obj;
		}
		
		return kNormalEval;
	}
}