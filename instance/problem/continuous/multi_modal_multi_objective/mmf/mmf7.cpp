#include "mmf7.h"
#include "../../../../../core/global.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include <fstream>

namespace ofec {
	void MMF7::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		m_domain.setRange(1., 3., 0);
		m_domain.setRange(-1., 1., 1);

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

	void MMF7::updateOptima() {
		m_optima.reset(new Optima<>);
		loadParetoFront(10000);
#ifdef OFEC_DEMO
		sampleParetoSets(10000);
#endif //
	}

	void MMF7::loadParetoFront(size_t sample_num) {
		std::stringstream os1;
		os1 << g_working_dir << "/instance/problem/continuous/multi_modal_multi_objective/mmf/data/" << m_name << "_obj.txt";
		auto ss = os1.str();
		std::ifstream infile(ss);
		if (!infile.is_open())
		{
			throw MyExcept("open PF file of MMF7 problem is fail");
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
				throw MyExcept("open PF file of MMF7 problem is fail");
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

	void MMF7::sampleParetoFront(size_t sample_num) {
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

	std::vector<Real> MMF7::createVar(const std::vector<Real>& s) const {
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> new_var;
		for (size_t i = 0; i < s.size(); ++i) {
			new_var.push_back(s[i]);
		}
		Real t = std::fabs(new_var[0]-2);
		Real pi = OFEC_PI;
		new_var[1] = (0.3*t*t*std::cos(24*pi*t+4*pi)+0.6*t)*std::sin(6 * pi * t + pi);

		return new_var;
	}

	void MMF7::sampleParetoSets(size_t sample_num) {
		size_t num = sample_num;
		std::vector<VariableVector<Real>> all_vars;
		
		Real pi = OFEC_PI;
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 0) {
					temp_var[j] = 1 + 2 * (Real)i / (num - 1.);//x1=1,3
				}
				else {
					Real t = std::fabs(temp_var[0] - 2);
					temp_var[j] = (0.3 * t * t * std::cos(24 * pi * t + 4 * pi) + 0.6 * t) * std::sin(6 * pi * t + pi);
				}
			}
			all_vars.emplace_back(temp_var);
		}
		for (size_t i = 0; i < all_vars.size(); ++i) {
			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
		}
		m_optima->setVariableGiven(true);
	}

	void MMF7::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Real t = std::fabs(x[0] - 2);
		Real pi = OFEC_PI;
		
		obj[0] = t;
		obj[1] = 1 - std::pow(t, 0.5) + std::pow(x[1] - (0.3*t*t*std::cos(24*pi*t + 4*pi)+0.6*t)*std::sin(6*pi*t+pi), 2);

	}
}