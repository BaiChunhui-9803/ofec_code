#include"f46_rastrigin_rastrigin.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../core/problem/continuous/continuous.h"

#include<iostream>
#include<fstream>
#include <string>

namespace ofec {
	void BiobjF46::initialize_()
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

		m_condition_number = 10.;
		m_beta = 0.2;
		loadRotation(sqrt(m_condition_number));

		//生成最优解集
		m_optima.reset(new Optima<>);
		compute_Xopt();
		//要在各函数的最优解之间进行网格化采样
		//generate_points();
		compute_Fopt();
		m_optima->setVariableGiven(true);
		m_optima->setObjectiveGiven(true);

	}


	void BiobjF46::evaluateObjective(Real* x, std::vector<Real>& obj)//x为指针首地址
	{
		//假设每个函数的最优解分别为（0,...）(1,...)
		obj.resize(m_number_objectives, 0);

		//评价f1-最优解为(-1,..)
		size_t i, j;
		Real tmp = 0., tmp2 = 0.;
		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmpvect[i] += m_rot[i][j] * (x[j] - dynamic_cast<Optima<>&>(*m_optima).variable(0)[j]);
			}
		}
		irregularize(tmpvect.data());
		asyemmetricalize(tmpvect.data(), m_beta);
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_linearTF[i][j] * tmpvect[j];
			}
		}
		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp += cos(2. * OFEC_PI * trasX[i]);
			tmp2 += trasX[i] * trasX[i];
		}
		obj[0] = 10. * ((Real)m_number_variables - tmp) + tmp2;

		//评价f2-最优解为(2,..)
		tmp = 0., tmp2 = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmpvect[i] += m_rot[i][j] * (x[j] - dynamic_cast<Optima<>&>(*m_optima).variable(1)[j]);
			}
		}
		irregularize(tmpvect.data());
		asyemmetricalize(tmpvect.data(), m_beta);
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_linearTF[i][j] * tmpvect[j];
			}
		}
		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp += cos(2. * OFEC_PI * trasX[i]);
			tmp2 += trasX[i] * trasX[i];
		}
		obj[1] = 10. * ((Real)m_number_variables - tmp) + tmp2;


	}

	void BiobjF46::compute_Xopt()
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
		Real x1_min = m_domain[0].limit.first;//dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0];
		Real x1_max = m_domain[0].limit.second; //dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[0] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[0];
		Real x2_min = m_domain[1].limit.first;//dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1];
		Real x2_max = m_domain[1].limit.second; //dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] > dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1] ? dynamic_cast<const Optima<>*>(m_optima.get())->variable(0)[1] : dynamic_cast<const Optima<>*>(m_optima.get())->variable(1)[1];
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


	void BiobjF46::compute_Fopt()
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

	void BiobjF46::irregularize(Real* x)
	{
		Real c1, c2, x_;
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (x[i] > 0) {
				c1 = 10;	c2 = 7.9;
			}
			else {
				c1 = 5.5;	c2 = 3.1;
			}
			if (x[i] != 0) {
				x_ = log(fabs(x[i]));
			}
			else x_ = 0;
			x[i] = sign(x[i]) * exp(x_ + 0.049 * (sin(c1 * x_) + sin(c2 * x_)));
		}
	}

	void BiobjF46::asyemmetricalize(Real* x, Real belta)
	{
		if (m_number_variables == 1) return;
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (x[i] > 0) {
				x[i] = pow(x[i], 1 + belta * i * sqrt(x[i]) / (m_number_variables - 1));
			}
		}
	}

	bool BiobjF46::loadRotation(Real base_)
	{
		computeRotation(m_rot, m_number_variables);
		computeRotation(m_rot2, m_number_variables);
		/* decouple scaling from function definition*/
		size_t i, j, k;
		m_linearTF.resize(m_number_variables, m_number_variables);
		for (i = 0; i < m_number_variables; ++i) {
			for (j = 0; j < m_number_variables; ++j) {
				m_linearTF[i][j] = 0.;
				for (k = 0; k < m_number_variables; ++k) {
					m_linearTF[i][j] += m_rot[i][k] * pow(base_, ((Real)k) / ((Real)(m_number_variables - 1))) * m_rot2[k][j];
				}
			}
		}
		return true;
	}

	void BiobjF46::reshape(Matrix& B, std::vector<Real>& vector, size_t m, size_t n)
	{
		B.resize(m, n);
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				B[i][j] = vector[j * m + i];
			}
		}
	}

	void BiobjF46::computeRotation(Matrix& rot, size_t Dim)
	{
		m_norRand.resize(Dim * Dim);
		for (auto& i : m_norRand) i = m_random.get()->normal.next();
		reshape(rot, m_norRand, Dim, Dim);
		/*1st coordinate is row, 2nd is column.*/
		size_t i, j, k;
		Real prod;
		for (i = 0; i < Dim; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				prod = 0;
				for (k = 0; k < Dim; ++k)
				{
					prod += rot[k][i] * rot[k][j];
				}
				for (k = 0; k < Dim; ++k)
				{
					rot[k][i] = rot[k][i] - prod * rot[k][j];
				}
			}
			prod = 0;
			for (k = 0; k < Dim; ++k)
			{
				prod += rot[k][i] * rot[k][i];
			}
			for (k = 0; k < Dim; ++k)
			{
				rot[k][i] = rot[k][i] / sqrt(prod);
			}
		}
	}

}