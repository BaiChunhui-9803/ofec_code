#include"f1_sphere_sphere.h"

#include<iostream>
#include<fstream>
#include <string>

namespace ofec {
	void BiobjF1::initialize_()
	{
		Continuous::initialize_();//调用problem：：initialize_,其中factory<problem>得到类似zdt1等具体函数
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
		compute_Fopt();
		m_optima->setVariableGiven(true);
		m_optima->setObjectiveGiven(true);

	}

	void BiobjF1::evaluateObjective(Real* x, std::vector<Real>& obj)//x为指针首地址
	{
		//假设每个函数的最优解分别为（0,...）(1,...)
		obj.resize(m_number_objectives, 0);

		for (size_t i = 0; i < m_number_objectives; i++)
		{
			for (size_t j = 0; j < m_number_variables; j++)
			{
				auto temp = dynamic_cast<const Optima<>*>(m_optima.get())->variable(i)[j];
				obj[i] += std::powf(x[j] - temp, 2);
			}
		}
	}

	void BiobjF1::compute_Xopt()
	{
		//目标1最优点
		VariableVector<Real> temp1(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			temp1[i] = 8 * floor(1e4 * m_random.get()->uniform.next()) / 1e4 - 4;
			if (temp1[i] == 0.0) {
				temp1[i] = -1e-5;
			}
		}
		dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp1);
		//目标2最优点
		VariableVector<Real> temp2(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			temp2[i] = 8 * floor(1e4 * m_random.get()->uniform.next()) / 1e4 - 4;
			if (temp2[i] == 0.0) {
				temp2[i] = -1e-5;
			}
		}
		dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp2);
		//把这两点间的ps录入到append var
		size_t sample_frequence = 1000;
		std::vector<std::vector<Real>> temp_sample(m_number_variables);
		Real min;
		Real max;
		for (size_t i = 0; i < m_number_variables; i++)
		{
			if (temp1[i] != temp2[i])
			{
				temp1[i] > temp2[i] ? (min = temp2[i], max = temp1[i]) : (min = temp1[i], max = temp2[i]);
				for (size_t j = 0; j < sample_frequence; j++)
				{
					temp_sample[i].push_back(m_random.get()->uniform.nextNonStd(min, max));
				}
				std::sort(temp_sample[i].begin(), temp_sample[i].end());//从小到大排序
			}
			else
			{
				for (size_t j = 0; j < sample_frequence; j++)
				{
					temp_sample[i].push_back(temp1[i]);
				}
			}
		}
		//插入最优解
		if ((std::signbit(temp1[0] - temp2[0]) == std::signbit(temp1[1] - temp2[1])))
		{
			for (size_t i = 0; i < sample_frequence; i++)
			{
				VariableVector<Real> temp(m_number_variables);
				temp[0] = temp_sample[0][i];
				temp[1] = temp_sample[1][i];
				dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp);
			}
		}
		else if ((std::signbit(temp1[0] - temp2[0]) != std::signbit(temp1[1] - temp2[1])))
		{
			std::sort(temp_sample[0].begin(), temp_sample[0].end(), std::greater<int>());
			for (size_t i = 0; i < sample_frequence; i++)
			{
				//dynamic_cast<Optima<>*>(m_optima.get())->appendVar({ temp_sample[0][i], temp_sample[1][i] });这种写法需要类型转换，我还不会//VariableVector<Real>
				VariableVector<Real> temp(m_number_variables);
				temp[0] = temp_sample[0][i];
				temp[1] = temp_sample[1][i];
				dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp);
			}
		}
		else if (dynamic_cast<Optima<>*>(m_optima.get())->numberVariables() < sample_frequence)
		{
			throw MyExcept("the number of optimal Variables may wrong.");//弹出警告
		}

	}

	void BiobjF1::compute_Fopt()
	{

		for (size_t i = 0; i < m_optima.get()->numberVariables(); i++)
		{
			std::vector<Real> temp(m_number_objectives);
			auto& a = (dynamic_cast<Optima<>*>(m_optima.get())->variable(i)[0]);
			evaluateObjective(&a, temp);
			m_optima->appendObj(temp);
		}
	}
}