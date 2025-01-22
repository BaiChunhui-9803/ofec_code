#include"f4_sphere_original_rosenbrock.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"

#include<iostream>
#include<fstream>
#include <string>

namespace ofec {
	void BiobjF4::initialize_()
	{
		Continuous::initialize_();//调用problem：：initialize_,其中factory<problem>得到zdt1等具体函数
		auto& v = *m_param;
		resizeObjective(v.get<int>("number of objectives"));
		if (m_number_objectives != 2)
			throw MyExcept("The number of objectives must be equal to 2");
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		m_optimize_mode[1] = OptimizeMode::kMinimize;

		resizeVariable(v.get<int>("number of variables"));//recomend n=2
		setDomain(-5, 5);
		//生成最优解集
		m_optima.reset(new Optima<>);
		compute_Xopt();
		//要在各函数的最优解之间进行网格化采样
		//generate_points();
		compute_Fopt();
		m_optima->setVariableGiven(true);
		m_optima->setObjectiveGiven(true);

	}


	void BiobjF4::evaluateObjective(Real* x, std::vector<Real>& obj)//x为指针首地址
	{
		//假设每个函数的最优解分别为（0,...）(1,...)
		obj.resize(m_number_objectives, 0);

		//评价f1-最优解为(-1,..)
		Real temp1 = -1.;
		for (size_t j = 0; j < m_number_variables; j++)
		{
			obj[0] += std::powf(x[j] - temp1, 2);
		}

		//评价f2-最优解为(2,..)
		Real temp2 = 2.;
		Real m_scales = 1.; Real f_opt = 0.; Real tmp;//copy from [paper] 2009 -real parameter.
		std::vector<Real> trasX(m_number_variables);
		for (size_t j = 0; j < m_number_variables; ++j) {
			trasX[j] = m_scales * (x[j] - temp2) + 1;
		}

		for (size_t j = 0; j < m_number_variables - 1; ++j)
		{
			tmp = (trasX[j] * trasX[j] - trasX[j + 1]);
			obj[1] += tmp * tmp;
		}
		obj[1] *= 1e2;
		for (size_t j = 0; j < m_number_variables - 1; ++j)
		{
			tmp = (trasX[j] - 1.);
			obj[1] += tmp * tmp;
		}

	}

	void BiobjF4::compute_Xopt()
	{
		//目标1最优点
		VariableVector<Real> temp1(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			temp1[i] = -1.;
		}
		dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp1);

		//目标2最优点
		VariableVector<Real> temp2(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			temp2[i] = 2.;
		}

		dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp2);

		//取各函数最优解之间评价
		Real x1_min = dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0];
		Real x1_max = dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0];
		Real x2_min = dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1];
		Real x2_max = dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1];
		//根据两个单目标的最优解产生一个box的sample
		for (size_t i = 0; i < sample_num_dimension; i++)
		{
			VariableVector<Real> temp(m_number_variables);
			Real x1 = x1_min + (Real)i / sample_num_dimension * (x1_max - x1_min);
			temp[0] = x1;
			for (size_t j = 0; j < sample_num_dimension; j++)
			{
				Real x2 = x2_min + (Real)j / sample_num_dimension * (x2_max - x2_min);
				temp[1] = x2;
				dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp);
			}
		}

	}


	void BiobjF4::compute_Fopt()
	{
		//trick把极值放在前两个，用来画boundary

		for (size_t i = 0; i < 2; i++)
		{
			std::vector<Real> temp(m_number_objectives);
			auto& a = (dynamic_cast<Optima<>*>(m_optima.get())->variable(i)[0]);
			evaluateObjective(&a, temp);
			m_optima->appendObj(temp);
		}

		if (m_number_objectives == 2)
		{
			std::vector<std::pair<Real, Real>> temp_optimal;
			std::pair<Real, Real>temp;
			temp.first = dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(0) < dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(0) ? dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(0) : dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(0);
			temp.second = dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(0) > dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(0) ? dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(0) : dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(0);
			temp_optimal.push_back(temp);

			std::pair<Real, Real>temp2;
			temp2.first = dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(1) < dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(1) ? dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(1) : dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(1);
			temp2.second = dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(1) > dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(1) ? dynamic_cast<Optima<>*>(m_optima.get())->objective(0).at(1) : dynamic_cast<Optima<>*>(m_optima.get())->objective(1).at(1);
			temp_optimal.push_back(temp2);

			m_optima->setObjRange(temp_optimal);
		}



		for (size_t i = 2; i < m_optima.get()->numberVariables(); i++)
		{
			std::vector<Real> temp(m_number_objectives);
			auto& a = (dynamic_cast<Optima<>*>(m_optima.get())->variable(i)[0]);
			evaluateObjective(&a, temp);
			m_optima->appendObj(temp);
		}

		//对这些解进行pareto 排序
		std::vector<std::vector<Real>> objs;
		for (size_t i = 0; i < m_optima->numberVariables(); i++)
			objs.emplace_back(m_optima->objective(i));
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);

		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] != 0)
			{
				dynamic_cast<Optima<>*>(m_optima.get())->eraseVar(i);
				dynamic_cast<Optima<>*>(m_optima.get())->eraseObj(i);
				rank.erase(rank.begin() + i);
				i--;
			}
		}
	}

}