#include "glt.h"
#include "../../../../../core/global.h"
#include <fstream>

namespace ofec {
	void GLT::initialize_() {
		Continuous::initialize_();
		auto &v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		m_num_reference_points = v.get<int>("# reference points");
		m_domain.setRange(0, 1, 0);
		if (m_name == "MOP_GLT5" || m_name == "MOP_GLT6") {
			m_domain.setRange(0, 1, 1);
			for (size_t i = 2; i < m_number_variables; i++)
				m_domain.setRange(-1, 1, i);
		}
		else {
			for (size_t i = 1; i < m_number_variables; i++)
				m_domain.setRange(-1, 1, i);
		}
		
		//m_domain_update = true;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOP_GLT5" || m_name == "MOP_GLT6") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
		}
		for (size_t i = 0; i < m_number_objectives; ++i) 
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoFront();		
	}

//	void GLT::updateOptima() {
//		m_optima.reset(new Optima<>);
//		loadParetoFront(10000);
//#ifdef OFEC_DEMO
//		sampleParetoSets(m_num_reference_points);
//#endif //
//	}
//
//	void GLT::loadParetoFront(size_t sample_num) {
//		std::stringstream os1;
//		os1 << g_working_dir << "/instance/problem/continuous/multi_objective/glt/data/" << m_name << "_obj.txt";
//		auto ss = os1.str();
//		std::ifstream infile(ss);
//		if (!infile.is_open())
//		{
//			throw MyExcept("open PF file of GLT problem is fail");
//		}
//		std::string str;
//		size_t line = 0;
//		while (getline(infile, str))
//			++line;
//		infile.close();
//		infile.clear();
//		if (line == 0) {
//			sampleParetoFront(10000);
//		}
//		else {
//			infile.open(os1.str());
//			if (!infile) {
//				throw MyExcept("open PF file of GLT problem is fail");
//			}
//			m_optima->resizeObjectiveSet(line);
//			for (size_t i = 0; i < line; i++) {
//				std::vector<Real> temp_obj(m_number_objectives);
//				for (size_t j = 0; j < m_number_objectives; j++)
//					infile >> temp_obj[j];
//				m_optima->setObjective(temp_obj, i);
//			}
//			m_optima->setObjectiveGiven(true);
//			infile.close();
//			infile.clear();
//		}
//	}
//
//	void GLT::sampleParetoFront(size_t sample_num) {
//		size_t num = sample_num;
//		std::vector<std::vector<Real>> all_objs;
//		if (m_name == "MOP_GLT1" || m_name == "MOP_GLT4") {
//			for (size_t i = 0; i < num; i++) {
//				std::vector<Real> temp_obj(m_number_objectives);
//				if (m_name == "MOP_GLT1") {
//					for (size_t j = 0; j < m_number_objectives; j++) {
//						if (j == 0)
//							temp_obj[j] = (Real)i / (num - 1.);
//						else
//							temp_obj[j] = 2 - temp_obj[0] - sign(std::cos(2 * OFEC_PI * temp_obj[0]));
//					}
//				}
//				else if (m_name == "MOP_GLT4") {
//					for (size_t j = 0; j < m_number_objectives; j++) {
//						if (j == 0)
//							temp_obj[j] = (Real)i / (num - 1.);
//						else
//							temp_obj[j] = 2 - 2 * std::pow(temp_obj[0],0.5) * std::pow(std::cos(3*OFEC_PI* temp_obj[0]* temp_obj[0]), 2);
//					}
//				}
//				all_objs.emplace_back(temp_obj);
//			}
//			//first nd-sort, then add into optima
//			std::vector<std::vector<Real>*> objs;
//			for (size_t i = 0; i < all_objs.size(); ++i) {
//				objs.emplace_back(&all_objs[i]);
//			}
//			std::vector<int> rank;
//			ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
//			for (size_t i = 0; i < rank.size(); ++i) {
//				if (rank[i] == 0)
//					m_optima->appendObj(all_objs[i]);
//			}
//		}
//		else if(m_name == "MOP_GLT2" || m_name == "MOP_GLT3"){
//			for (size_t i = 0; i < num; i++) {
//				std::vector<Real> temp_obj(m_number_objectives);
//				Real temp = (Real)i / (num - 1.);
//				for (size_t j = 0; j < m_number_objectives; j++) {
//					if (m_name == "MOP_GLT2") {
//						if (j == 0) 
//							temp_obj[j] = 1 - std::cos(temp * OFEC_PI / 2);
//						else
//							temp_obj[j] = 10 - 10*std::sin(temp * OFEC_PI / 2);
//					}
//					else {
//						if (j == 0) 
//							temp_obj[j] = temp;
//						else {
//							if (temp<=0.05)
//								temp_obj[j] = 1 - 19*temp_obj[0];
//							else 
//								temp_obj[j] = 1./19 - temp_obj[0]/19.;
//						}
//					}
//				}
//				all_objs.emplace_back(temp_obj);
//			}
//			for (size_t i = 0; i < all_objs.size(); ++i) {
//				m_optima->appendObj(all_objs[i]);
//			}
//		}
//		else {
//			size_t div_num = 100;
//			for (size_t i = 0; i < div_num; i++) {
//				Real temp1= (Real)i / (div_num - 1.);
//				for (size_t j = 0; j < div_num; j++) {
//					Real temp2= (Real)j / (div_num - 1.);
//					std::vector<Real> temp_obj(m_number_objectives);
//					temp_obj[0] = (1-std::cos(temp1*OFEC_PI/2))*(1 - std::cos(temp2 * OFEC_PI / 2));
//					temp_obj[1] = (1 - std::cos(temp1 * OFEC_PI / 2)) * (1 - std::sin(temp2 * OFEC_PI / 2));
//					if (m_name == "MOP_GLT5") 
//						temp_obj[2] = 1 - std::sin(temp1 * OFEC_PI / 2);
//					else 
//						temp_obj[2] = 2 - std::sin(temp1 * OFEC_PI / 2)-sign(std::cos(4*temp1*OFEC_PI));
//					all_objs.emplace_back(temp_obj);
//				}
//			}
//			if (m_name == "MOP_GLT6") {
//				//先排序，再加入
//				std::vector<std::vector<Real>*> objs;
//				for (size_t i = 0; i < all_objs.size(); ++i) {
//					objs.emplace_back(&all_objs[i]);
//				}
//				std::vector<int> rank;
//				ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
//				for (size_t i = 0; i < rank.size(); ++i) {
//					if (rank[i] == 0) {
//						m_optima->appendObj(all_objs[i]);
//					}
//				}
//			}
//			else {
//				for (size_t i = 0; i < all_objs.size(); ++i) {
//					m_optima->appendObj(all_objs[i]);
//				}
//			}
//		}
//		m_optima->setObjectiveGiven(true);
//
//		//write into file
//		std::ofstream out1;
//		std::stringstream os1;
//		os1 << g_working_dir << "/instance/problem/continuous/multi_objective/glt/data/" << m_name << "_obj.txt";
//		out1.open(os1.str());
//		if (!out1) {
//			throw("can not open file");
//		}
//		for (size_t i = 0; i < m_optima->numberObjectives(); ++i) {
//			auto& objs = m_optima->objective(i);
//			for (size_t j = 0; j < objs.size(); ++j) {
//				out1 << objs[j] << " ";
//			}
//			out1 << std::endl;
//		}
//		out1.close();
//	}

	std::vector<Real> GLT::createVar(const std::vector<Real>& s) const {
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> new_var;
		for (size_t i = 0; i < s.size(); ++i) {
			new_var.push_back(s[i]);
		}

		if (m_name != "GLT6" && m_name != "GLT5") {
			
			if (m_name == "GLT1" || m_name == "GLT4") {
				
			}
			else {
				for (size_t j = 1; j < m_number_variables; j++) {
					new_var[j] = std::sin(2 * OFEC_PI * new_var[0] + (1 + j) * OFEC_PI / m_number_variables);
				}
			}
		}
		else if(m_name == "GLT5"){
			for (size_t k = 2; k < m_number_variables; k++) {
				new_var[k] = 2 * new_var[1] * std::sin(2 * OFEC_PI * new_var[0] + (1 + k) * OFEC_PI / m_number_variables);
			}
		}
		else {
			
		}
		return new_var;
	}
}