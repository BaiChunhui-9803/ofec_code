#include "uf.h"
#include <fstream>
#include "../../../../../core/global.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"

namespace ofec {
	void UF::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;;
		resizeVariable(v.get<int>("number of variables"));  //recomend n=10
		if (m_name == "MOP_UF03")
			setDomain(0., 1.);
		else {
			if (m_name == "MOP_UF04") {
				m_domain.setRange(0., 1., 0);
				for (size_t i = 1; i < m_number_variables; ++i)
					m_domain.setRange(-2., 2., i);
			}
			else if (m_name == "MOP_UF08" || m_name == "MOP_UF09" || m_name == "MOP_UF10") {
				for (size_t i = 0; i < m_number_variables; ++i) {
					if (i == 0 || i == 1)
						m_domain.setRange(0., 1., i);
					else
						m_domain.setRange(-2., 2., i);
				}
			}
			else {
				m_domain.setRange(0., 1., 0);
				for (size_t i = 1; i < m_number_variables; ++i)
					m_domain.setRange(-1., 1., i);
			}
			m_domain_update = true;
		}

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOP_UF08" || m_name == "MOP_UF09" || m_name == "MOP_UF10") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
			if (m_number_variables < 3)
				throw MyExcept("The number of variables must be no less than 3.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
			/*if (m_name == "MOP_UF03" && m_number_variables < 3)
				throw MyExcept("The number of variables must be no less than 3.");*/
		}
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoFront();
	}

	void UF::updateOptima() {
		m_optima.reset(new Optima<>);
		loadParetoFront(10000);
		sampleParetoSets(10000);
	}

	void UF::loadParetoFront(size_t sample_num) {
		std::stringstream os1;
		os1 << g_working_dir << "/instance/problem/continuous/multi_objective/uf/data/" << m_name << "_obj.txt";
		auto ss = os1.str();
		std::ifstream infile(ss);
		if (!infile.is_open())
		{
			throw MyExcept("open PF file of UF problem is fail");
		}
		std::string str;
		size_t line = 0;
		while (getline(infile, str))
			++line;
		infile.close();
		infile.clear();
		if (line == 0) {
			sampleParetoFront(10000);
		}
		else {
			infile.open(os1.str());
			if (!infile) {
				throw MyExcept("open PF file of UF problem is fail");
			}
			m_optima->resizeObjectiveSet(line);
			for (size_t i = 0; i < line; i++) {
				std::vector<Real> temp_obj(m_number_objectives);
				for (size_t j = 0; j < m_number_objectives; j++)
					infile >> temp_obj[j];
				m_optima->setObjective(temp_obj, i);
			}
			m_optima->setObjectiveGiven(true);
			infile.close();
			infile.clear();
		}
	}

	void UF::sampleParetoFront(size_t sample_num) {
		size_t num = sample_num;
		std::vector<std::vector<Real>> all_objs;
		if (m_name == "MOP_UF01" || m_name == "MOP_UF02" || m_name == "MOP_UF03") {
			for (size_t i = 0; i < num; i++) {
				std::vector<Real> temp_obj(m_number_objectives);
				temp_obj[0] = (Real)i / (num - 1.);
				temp_obj[1] = 1 - std::sqrt(temp_obj[0]);
				all_objs.emplace_back(temp_obj);
			}
		}
		else if (m_name == "MOP_UF04") {
			for (size_t i = 0; i < num; i++) {
				std::vector<Real> temp_obj(m_number_objectives);
				temp_obj[0] = (Real)i / (num - 1.);
				temp_obj[1] = 1 - temp_obj[0] * temp_obj[0];
				all_objs.emplace_back(temp_obj);
			}
		}
		else if (m_name == "MOP_UF05") {
			size_t N = 10;
			for (size_t i = 0; i < 2*N; i++) {
				std::vector<Real> temp_obj(m_number_objectives);
				temp_obj[0] = (Real)i / 2./N;
				temp_obj[1] = 1. - (Real)i / 2. / N;
				all_objs.emplace_back(temp_obj);
			}
		}
		else if (m_name == "MOP_UF06") {
			size_t N = 2;
			for (size_t i = 0; i < N; i++) {
				size_t samples = 2000;
				Real start_p = (2 * (i+1) - 1) / 2. / N;
				Real end_p = (i+1.) / N;
				for (size_t j = 0; j < samples; ++j) {
					std::vector<Real> temp_obj(m_number_objectives);
					temp_obj[0] = start_p + j * (end_p-start_p) / (samples-1);
					temp_obj[1] = 1. - temp_obj[0];
					all_objs.emplace_back(temp_obj);
				}
			}
			std::vector<Real> temp{ 0.,1. };
			all_objs.emplace_back(temp);
		}
		else if (m_name == "MOP_UF07") {
			for (size_t i = 0; i < num; i++) {
				std::vector<Real> temp_obj(m_number_objectives);
				temp_obj[0] = (Real)i / (num - 1.);
				temp_obj[1] = 1 - temp_obj[0];
				all_objs.emplace_back(temp_obj);
			}
		}
		else if (m_name == "MOP_UF08"|| m_name == "MOP_UF10") {
			size_t div_num = 100;
			for (size_t i = 0; i < div_num; i++) {
				Real temp1 = (Real)i / (div_num - 1.);
				for (size_t j = 0; j < div_num; j++) {
					Real temp2 = (Real)j / (div_num - 1.);
					std::vector<Real> temp_obj(m_number_objectives);
					temp_obj[0] = std::cos(temp1 * OFEC_PI / 2) * std::cos(temp2 * OFEC_PI / 2);
					temp_obj[1] = std::cos(temp1 * OFEC_PI / 2) * std::sin(temp2 * OFEC_PI / 2);
					temp_obj[2] = std::sin(temp1 * OFEC_PI / 2);
					all_objs.emplace_back(temp_obj);
				}
			}
		}
		else if (m_name == "MOP_UF09") {
			size_t div_num = 100;
			Real e = 0.1;
			std::vector<std::vector<Real>> temp_all_objs;
			for (size_t i = 0; i < div_num; i++) {
				Real temp1 = (Real)i / (div_num - 1.);
				for (size_t j = 0; j < div_num; j++) {
					Real temp2 = (Real)j / (div_num - 1.);
					std::vector<Real> temp_obj(m_number_objectives);
					Real temp3 = (1 + e) * (1 - 4 * std::pow(2 * temp1 - 1, 2));
					Real temp4 = temp3 > 0 ? temp3 : 0;
					temp_obj[0] = 0.5 *(2*temp1+temp4)*temp2;
					temp_obj[1] = 0.5 * (2+temp4-2*temp1) * temp2;
					temp_obj[2] = 1-temp2;
					temp_all_objs.emplace_back(temp_obj);
				}
			}
			// nd-sort
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < temp_all_objs.size(); ++i) {
				objs.emplace_back(&temp_all_objs[i]);
			}
			std::vector<int> rank;
			ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
			for (size_t i = 0; i < rank.size(); ++i) {
				if (rank[i] == 0) {
					all_objs.emplace_back(temp_all_objs[i]);
				}
			}
		}
		for (size_t i = 0; i < all_objs.size(); ++i) {
			m_optima->appendObj(all_objs[i]);
		}
		m_optima->setObjectiveGiven(true);

		//write into file
		std::ofstream out1;
		std::stringstream os1;
		os1 << g_working_dir << "/instance/problem/continuous/multi_objective/uf/data/" << m_name << "_obj.txt";
		out1.open(os1.str());
		if (!out1) {
			throw("can not open file");
		}
		for (size_t i = 0; i < m_optima->numberObjectives(); ++i) {
			auto& objs = m_optima->objective(i);
			for (size_t j = 0; j < objs.size(); ++j) {
				out1 << objs[j] << " ";
			}
			out1 << std::endl;
		}
		out1.close();
	}

	void UF::sampleParetoSets(size_t sample_num) {
		size_t num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		if (m_name == "MOP_UF01" || m_name == "MOP_UF04" || m_name == "MOP_UF07") {
			for (size_t i = 0; i < num; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				temp_var[0] = (Real)i / (num - 1.);
				for (size_t j = 1; j < m_number_variables; ++j) {
					temp_var[j] = std::sin(6 * OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				all_vars.emplace_back(temp_var);
			}
		}
		else if (m_name == "MOP_UF02") {
			for (size_t i = 0; i < num; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				temp_var[0] = (Real)i / (num - 1.);
				for (size_t j = 1; j < m_number_variables; ++j) {
					if(j%2==1)
						temp_var[j] = (0.3*std::pow(temp_var[0],2)*std::cos(24 * OFEC_PI * temp_var[0] + 4*(j + 1) * OFEC_PI / m_number_variables)+0.6*temp_var[0]) * std::sin(6* OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables);
					else
						temp_var[j] = (0.3 * std::pow(temp_var[0], 2) * std::cos(24 * OFEC_PI * temp_var[0] + 4 * (j + 1) * OFEC_PI / m_number_variables) + 0.6 * temp_var[0]) * std::cos(6 * OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				all_vars.emplace_back(temp_var);
			}
		}
		else if (m_name == "MOP_UF03") {
			for (size_t i = 0; i < num; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				temp_var[0] = (Real)i / (num - 1.);
				for (size_t j = 1; j < m_number_variables; ++j) {
					if(m_number_variables==2)
						temp_var[j] = std::pow(temp_var[0], 0.5 * (1.+3.*j/(m_number_variables-1.)));
					else
						temp_var[j] = std::pow(temp_var[0], 0.5 * (1. + 3. * (j - 1.) / (m_number_variables - 2.)));
				}
				all_vars.emplace_back(temp_var);
			}
		}
		else if (m_name == "MOP_UF05") {
			size_t N = 10;
			for (size_t i = 0; i < 2 * N; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				temp_var[0] = (Real)i / 2. / N;
				for (size_t j = 1; j < m_number_variables; ++j) {
					temp_var[j] = std::sin(6*OFEC_PI*temp_var[0]+(j+1)*OFEC_PI/m_number_variables);
				}
				all_vars.emplace_back(temp_var);
			}
		}
		else if (m_name == "MOP_UF06") {
			std::vector<VariableVector<Real>> temp_all_vars;
			for (size_t i = 0; i < num; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				temp_var[0] = (Real)i / (num - 1.);
				for (size_t j = 1; j < m_number_variables; ++j) {
					temp_var[j] = std::sin(6 * OFEC_PI * temp_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				temp_all_vars.emplace_back(temp_var);
			}
			//nd-sort
			for (size_t i = 0; i < temp_all_vars.size(); ++i) {
				std::vector<Real> temp_obj(m_number_objectives);
				evaluateObjective(temp_all_vars[i].data(), temp_obj);
				all_objs.emplace_back(temp_obj);
			}
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < all_objs.size(); ++i) {
				objs.emplace_back(&all_objs[i]);
			}
			std::vector<int> rank;
			ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
			for (size_t i = 0; i < rank.size(); ++i) {
				if (rank[i] == 0) {
					all_vars.emplace_back(temp_all_vars[i]);
				}
			}
		}
		else if (m_name == "MOP_UF08" || m_name == "MOP_UF10") {
			size_t div_num = 100;
			for (size_t i = 0; i < div_num; i++) {
				Real temp1 = (Real)i / (div_num - 1.);
				for (size_t j = 0; j < div_num; j++) {
					Real temp2 = (Real)j / (div_num - 1.);
					VariableVector<Real> temp_var(m_number_variables);
					temp_var[0] = temp1;
					temp_var[1] = temp2;
					for (size_t k = 2; k < m_number_variables; ++k) {
						temp_var[k] = 2 * temp2 * std::sin(2*OFEC_PI*temp1+(k+1)*OFEC_PI/m_number_variables);
					}
					all_vars.emplace_back(temp_var);
				}
			}
		}
		else if (m_name == "MOP_UF09") {
			size_t div_num = 100;
			for (size_t i = 0; i < div_num; i++) {
				Real temp1 = (Real)i / (div_num - 1.);
				if ((temp1 >= 0. && temp1 <= 0.25) || (temp1 >= 0.75 && temp1 <= 1.)) {
					for (size_t j = 0; j < div_num; j++) {
						Real temp2 = (Real)j / (div_num - 1.);
						VariableVector<Real> temp_var(m_number_variables);
						temp_var[0] = temp1;
						temp_var[1] = temp2;
						for (size_t k = 2; k < m_number_variables; ++k) {
							temp_var[k] = 2 * temp2 * std::sin(2 * OFEC_PI * temp1 + (k + 1) * OFEC_PI / m_number_variables);
						}
						all_vars.emplace_back(temp_var);
					}
				}
			}
		}
		for (size_t i = 0; i < all_vars.size(); ++i) {
			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
		}
		m_optima->setVariableGiven(true);
	}

	//void UF::loadParetoFront() {
	//	std::stringstream os;
	//	os << g_working_dir << "/instance/problem/continuous/multi_objective/uf/data/" << m_name << ".dat";
	//	std::ifstream infile(os.str());
	//	if (!infile) {
	//		//std::cout << "open file is failed" << std::endl;
	//		throw "open PF file of UF problem is fail";
	//		return;
	//	}
	//	std::string str;
	//	size_t line = 0;
	//	while (getline(infile, str))
	//		++line;
	//	m_optima.reset(new Optima<>());
	//	m_optima->resizeObjectiveSet(line);
	//	infile.close();
	//	infile.clear();
	//	infile.open(os.str());
	//	for (size_t i = 0; i < line; i++) {
	//		std::vector<Real> temp_obj(m_number_objectives);
	//		for (size_t j = 0; j < m_number_objectives; j++)
	//			infile >> temp_obj[j];
	//		m_optima->setObjective(temp_obj, i);
	//	}
	//	m_optima->setObjectiveGiven(true);
	//	infile.close();
	//}
}