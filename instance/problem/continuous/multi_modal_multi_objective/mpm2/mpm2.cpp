#include"mpm2.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"

namespace ofec {

	void MPM2::initialize_()
	{
		Continuous::initialize_();//调用problem：：initialize_,其中factory<problem>得到类似zdt1等具体函数
		auto& v =*m_param;


		mpm_objs.emplace_back(new MultiplePeaksModel(v.get<int>("number of variables"),
			v.get<int>("obj-1.n.peaks"), v.get<int>("obj-1.n.optima"), v.get<std::string>("obj-1.topology"),
			v.get<int>("obj-1.correlation")));

		mpm_objs.emplace_back(new MultiplePeaksModel(v.get<int>("number of variables"),
			v.get<int>("obj-2.n.peaks"), v.get<int>("obj-2.n.optima"), v.get<std::string>("obj-2.topology"),
			v.get<int>("obj-2.correlation")));//var,peaksum,optima,topology,correlation



		resizeObjective(v.get<int>("number of objectives"));
		if (m_number_objectives != 2)
			throw MyExcept("The objectives of this version MPM2  of objectives set 2, but it can be 3+, i did not write");
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		m_optimize_mode[1] = OptimizeMode::kMinimize;

		resizeVariable(v.get<int>("number of variables"));//recomend n=2
		setDomain(0, 1);
		//生成最优解集
		m_optima.reset(new Optima<>);
		//--为了moPLOT的可视化，这里先将下述求最优解的step屏蔽掉
		compute_opt();
		m_optima->setVariableGiven(true);
		m_optima->setObjectiveGiven(true);
		
	}

	void MPM2::evaluateObjective(Real* x, std::vector<Real>& obj)
	{
		for (size_t i = 0; i < mpm_objs.size(); i++)
		{
			obj[i] = 1.0 - mpm_objs[i].get()->getGvalue(x);
		}

	}

	void MPM2::compute_opt()
	{
		size_t sample_frequence = 100;
		Real x1_min = m_domain[0].limit.first;
		Real x1_max = m_domain[0].limit.second;
		Real x2_min = m_domain[1].limit.first;
		Real x2_max = m_domain[1].limit.second;
		//采样
		for (size_t i = 0; i < sample_frequence; i++)
		{
			VariableVector<Real> temp(m_number_variables);

			Real x1 = x1_min + (Real)i / sample_frequence * (x1_max - x1_min);
			temp[0] = x1;
			for (size_t j = 0; j < sample_frequence; j++)
			{
				Real x2 = x2_min + (Real)j / sample_frequence * (x2_max - x2_min);
				temp[1] = x2;
				dynamic_cast<Optima<>*>(m_optima.get())->appendVar(temp);
			}
		}
		//评价
		for (size_t i = 0; i < m_optima.get()->numberVariables(); i++)
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