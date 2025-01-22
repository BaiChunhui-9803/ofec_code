#include "MPM.h"

namespace ofec {

	MultiplePeaksModel::MultiplePeaksModel(const int vars, const int peaksnum, const int optimanum, const std::string topology, const int correlation, const std::string peakShape, const bool rotation)//:m_peaks(peaksnum)
	{
		if (optimanum <= 0)
		{
			throw"optima number must bigger than 0";
		}
		num_var = vars;
		num_peaks = peaksnum;
		num_optimal = optimanum;
		m_topology = topology;
		shapeHeightCorrelation = correlation;
		m_peakShape = peakShape;
		m_rotation = rotation;

		gen = std::mt19937(rd());
		uni_domain = std::uniform_real_distribution<double>(M_domain_min, M_domain_max);
		nor_domain = std::normal_distribution<double>(M_domain_min, M_domain_max);
		uni_rand_shapeRange = std::uniform_real_distribution<double>(1.5, 2.5);
		uni_rand_radiusRange = std::uniform_real_distribution<double>(0.25 * std::sqrt(vars), 0.5 * std::sqrt(vars));
		uni_rand_heightRange = std::uniform_real_distribution<double>(0.5, 0.99);

		for (int i = 0; i < num_optimal; i++) {
			double temp_shape = uni_rand_shapeRange(gen);
			double temp_radius = uni_rand_radiusRange(gen);

			m_peaks.emplace_back(new Peak(num_var, 1.0, temp_shape, temp_radius, {}, m_rotation, m_peakShape));

		}
		createInstanceWithExactNumberOfOptima();
	}


	void MultiplePeaksModel::randomUniformPeaks(int num_peaks)
	{

		if (num_peaks > 0)
		{
			for (int i = 0; i < num_peaks; i++) {
				double temp = uni_rand_heightRange(gen);
				while (temp == M_height) { temp = uni_rand_heightRange(gen); }
				m_peaks.push_back(std::make_unique<Peak>(num_var, temp, uni_rand_shapeRange(gen), uni_rand_radiusRange(gen)));
			}
		}
	}


	void MultiplePeaksModel::createInstance()
	{
		if (m_topology == "funnel")
		{
			if (num_optimal == 1)
			{
				makeFunnel();
			}

		}
		if (shapeHeightCorrelation != 0)
		{
			//按高度排
			std::sort(m_peaks.begin(), m_peaks.end(), [](const std::unique_ptr<Peak>& p1, const std::unique_ptr<Peak>& p2) {return p1->m_height > p2->m_height; });
			std::vector<double> temp_shape;
			for (size_t i = 0; i < m_peaks.size(); i++)
			{
				temp_shape.push_back(m_peaks.at(i).get()->m_shape);
			}

			if (shapeHeightCorrelation > 0)
			{
				std::sort(temp_shape.begin(), temp_shape.end(), std::greater<double>());
			}
			else
			{
				std::sort(temp_shape.begin(), temp_shape.end(), std::less<double>());
			}

			//按correlation 排  lambda表达式允许捕获局部变量，但是数据成员不是局部变量。可以捕获this使用类的数据成员。为+-由高到低
			for (size_t i = 0; i < m_peaks.size(); i++)
			{
				m_peaks[i].get()->m_shape = temp_shape[i];
			}
		}

	}

	void MultiplePeaksModel::makeFunnel()
	{
		// 排高度 高-低
		std::sort(m_peaks.begin(), m_peaks.end(), [](const std::unique_ptr<Peak>& p1, const std::unique_ptr<Peak>& p2) {return p1->m_height > p2->m_height; });
		std::vector<double> temp_height;
		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			temp_height.push_back(m_peaks[i].get()->m_height);
			m_peaks[i].get()->m_distance = 0.0;
		}


		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			double temp = 0.0;
			for (size_t j = 0; j < num_var; j++)
			{
				temp += std::pow(m_peaks[i]->m_position[j] - m_peaks[0]->m_position[j], 2);
			}
			m_peaks[i]->m_distance = temp;
		}
		// 排距离  小-大
		std::sort(m_peaks.begin(), m_peaks.end(), [](const std::unique_ptr<Peak>& p1, const std::unique_ptr<Peak>& p2) {return p1->m_distance < p2->m_distance; });
		//高度 高-低
		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			m_peaks[i].get()->m_height = temp_height[i];
		}

	}

	void MultiplePeaksModel::clusteredPeaks(int num_peaks)
	{
		//随机产生一个聚类中心
		std::vector<double> clusterCenter = m_peaks.at(0).get()->m_position;

		//围绕聚类中心，均匀分布峰
		for (int i = 0; i < num_peaks; i++) {
			std::vector<double> position(num_var);
			for (int j = 0; j < num_var; j++) {
				position[j] = nor_domain(gen) / 6.0 * std::sqrt(num_var) + clusterCenter[j];
				while (position[j] < M_domain_min || M_domain_max < position[j]) {
					if (M_domain_max < position[j]) {
						position[j] = M_domain_max - (position[j] - M_domain_max);
					}
					else if (position[j] < M_domain_min) {
						position[j] = -position[j];
					}
				}
			}
			m_peaks.push_back(std::make_unique<Peak>(num_var, uni_rand_heightRange(gen), uni_rand_shapeRange(gen), uni_rand_radiusRange(gen), position));
		}

	}



	double MultiplePeaksModel::G(Real* x, const std::unique_ptr<Peak>& peak)
	{
		//Mahalanobis distance
		std::vector<double> temp;
		for (size_t i = 0; i < num_var; i++)
		{
			temp.push_back(x[i] - peak.get()->m_position.at(i));
		}
		Eigen::VectorXd differenceVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(temp.data(), temp.size());

		double Mah_dis = std::sqrt(differenceVector.transpose() * peak.get()->m_covariance.inverse() * differenceVector);
		//求g
		return peak.get()->m_height / (1.0 + std::pow(Mah_dis, peak.get()->m_shape) / peak.get()->m_radius);
	}

	int MultiplePeaksModel::getActivePeak(Real* x)
	{
		//衡量该点与所有峰的目标函数值
		std::vector<double> count;
		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			count.push_back(G(x, m_peaks[i]));
		}
		return std::max_element(count.begin(), count.end()) - count.begin();
	}

	double MultiplePeaksModel::getGvalue(Real* x)
	{
		std::vector<double> count;
		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			count.push_back(G(x, m_peaks[i]));
		}
		return *std::max_element(count.begin(), count.end());
	}

	int MultiplePeaksModel::getLocalOptima()
	{
		int num_localOptima = 0;
		for (size_t i = 0; i < m_peaks.size(); i++)
		{
			if (getActivePeak(&m_peaks[i].get()->m_position[0]) == i) {
				num_localOptima++;
			}
			else
			{
				m_peaks.erase(m_peaks.begin() + i);
				m_peaks.shrink_to_fit();//for (auto& p : m_peaks)按照capacity
				i--;
			}


		}
		return num_localOptima;
	}

	void MultiplePeaksModel::createInstanceWithExactNumberOfOptima()
	{
		if (num_optimal > num_peaks)
		{
			throw "the number of peaks should larger than optima";
			return;
		}
		if (num_optimal == num_peaks)
		{
			//randomUniformPeaks(num_optimal-1);
			//createInstance();
			return;
		}
		//处理其他峰
		int remain_optimal = num_peaks - num_optimal;
		if (m_topology == "random")
		{
			randomUniformPeaks(remain_optimal);
		}
		else if (m_topology == "funnel")
		{
			clusteredPeaks(remain_optimal);
		}
		createInstance();
		//判断已生成的局部最优数目--是否与想要的一致。
		//通过采样，判断获得有效峰值的个数是否==获得局部最优解的数量,
		//即判断峰有没有被遮盖，即，该peak的position处的点，最大值是否为自己
		int current_num_optimal = getLocalOptima();
		int demand_optima = num_peaks - current_num_optimal;


		while (current_num_optimal < num_peaks * M_factor)
		{
			for (size_t i = 0; i < m_peaks.size(); i++)
			{
				m_peaks[i].get()->m_radius = m_peaks[i].get()->m_radius * M_radii;
			}
			//for (auto& p : m_peaks)按照capacity
			if (m_topology == "random")
			{
				randomUniformPeaks(demand_optima);
			}
			else if (m_topology == "funnel")
			{
				clusteredPeaks(demand_optima);
			}
			current_num_optimal = getLocalOptima();
		}


		while (demand_optima > 0)
		{
			if (m_topology == "random")
			{
				randomUniformPeaks(1);
			}
			else if (m_topology == "funnel")
			{
				clusteredPeaks(1);
			}
			demand_optima = num_peaks - getLocalOptima();
		}


		if (demand_optima <= 0)
		{
			createInstance();
		}

	}

}