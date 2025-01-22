#include "mmf3.h"
#include "../../../../../core/global.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include <fstream>

namespace ofec {
	void MMF3::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		m_domain.setRange(0., 1., 0);
		m_domain.setRange(0., 1.5, 1);

		m_domain_update = true;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_number_objectives != 2)
			throw MyExcept("The number of objectives must be equal to 2.");
		if (m_number_variables != 2)
			throw MyExcept("The number of variables must be equal to 2.");
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoFront();		
	}

	void MMF3::updateOptima() {
		m_optima.reset(new Optima<>);
		loadParetoFront(10000);
#ifdef OFEC_DEMO
		sampleParetoSets(10000);
#endif //
	}

	void MMF3::loadParetoFront(size_t sample_num) {
		std::stringstream os1;
		os1 << g_working_dir << "/instance/problem/continuous/multi_modal_multi_objective/mmf/data/" << m_name << "_obj.txt";
		auto ss = os1.str();
		std::ifstream infile(ss);
		if (!infile.is_open())
		{
			throw MyExcept("open PF file of MMF3 problem is fail");
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
				throw MyExcept("open PF file of MMF3 problem is fail");
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

	void MMF3::sampleParetoFront(size_t sample_num) {
		size_t num = sample_num;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < num; i++) {
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_objectives; j++) {
				if (j == 0)
					temp_obj[j] = (Real)i / (num - 1.);
				else
					temp_obj[j] = 1 - std::pow(temp_obj[0], 0.5);
			}
			all_objs.emplace_back(temp_obj);
		}
		for (size_t i = 0; i < all_objs.size(); ++i) {
			m_optima->appendObj(all_objs[i]);
		}
		m_optima->setObjectiveGiven(true);

		//write into file
		std::ofstream out1;
		std::stringstream os1;
		os1 << g_working_dir << "/instance/problem/continuous/multi_modal_multi_objective/mmf/data/" << m_name << "_obj.txt";
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

	std::vector<Real> MMF3::createVar(const std::vector<Real>& s) const {
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> new_var;
		for (size_t i = 0; i < s.size(); ++i) {
			new_var.push_back(s[i]);
		}
		if ((new_var[1] >= 0 && new_var[1] <= 0.5) ||( (new_var[1]>0.5 && new_var[1]<1) && (new_var[0] > 0.25 && new_var[0] < 1))) {
			new_var[0] = std::pow(new_var[1], 2);
		}
		else {
			new_var[0] = std::pow(new_var[1] - 0.5, 2);
		}

		return new_var;
	}

	void MMF3::sampleParetoSets(size_t sample_num) {
		size_t num = sample_num/2;
		std::vector<VariableVector<Real>> all_vars;
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 1) {
					temp_var[j] = (Real)i / (num - 1.);//x2=0-1
				}
				else {
					temp_var[j] = std::pow(temp_var[1], 2);
				}
			}
			all_vars.emplace_back(temp_var);
		}
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 1) {
					temp_var[j] = 0.5+(Real)i / (num - 1.);
				}
				else {
					temp_var[j] = std::pow(temp_var[1] - 0.5, 2);
				}
			}
			all_vars.emplace_back(temp_var);
		}
		for (size_t i = 0; i < all_vars.size(); ++i) {
			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
		}
		m_optima->setVariableGiven(true);
	}

	void MMF3::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real g = std::fabs(x[0] - 2.);

		obj[0] = x[0];
		if ((x[1] >= 0 && x[1] <= 0.5) || ((x[1] > 0.5 && x[1] < 1) && (x[0] > 0.25 && x[0] < 1))) {
			obj[1] = 1 - std::sqrt(x[0]) + 2 * (4 * std::pow(x[1] - std::sqrt(x[0]), 2) - 2 * std::cos(20 * (x[1] - std::sqrt(x[0])) * OFEC_PI / std::sqrt(2)) + 2);
		}
		else {
			obj[1] = 1 - std::sqrt(x[0]) + 2 * (4 * std::pow(x[1] - 0.5 - std::sqrt(x[0]), 2) - std::cos(20 * (x[1] - 0.5 - std::sqrt(x[0])) * OFEC_PI / std::sqrt(2)) + 2);
		}
		
	}
}