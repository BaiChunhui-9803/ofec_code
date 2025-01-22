#include "dpc.h"

namespace ofec {
	void DPC::getdc() {
		//double max_dis, min_dis;
		//max_dis = min_dis = mv_dis[0][1];
		//for (decltype(mv_dis.size()) i = 0; i < mv_dis.size() - 1; i++) {
		//	for (decltype(mv_dis.size()) j = i + 1; j < mv_dis.size(); j++) {
		//		if (mv_dis[i][j] > max_dis)
		//			max_dis = mv_dis[i][j];
		//		if (mv_dis[i][j] < min_dis)
		//			min_dis = mv_dis[i][j];
		//	}
		//}
		//m_dc = min_dis + (0.2 * max_dis);
		double percent = 10.0;
		int n = mv_dis.size() * (mv_dis.size() - 1) / 2;
		int position = static_cast<int>(n * percent / 100. + 0.5);
		std::multiset<double> sda;
		for (decltype(mv_dis.size()) i = 0; i < mv_dis.size() - 1; i++)
			for (decltype(mv_dis.size()) j = i + 1; j < mv_dis.size(); j++)
				sda.insert(mv_dis[i][j]);
		std::multiset<double>::iterator iter = sda.begin();
		n = 0;
		while (n < position - 1) {
			++iter;
			++n;
		}
		m_dc = *iter;
	}

	void DPC::getLocalDensity(LocalDensity LocalDensityVersion) {
		if (LocalDensityVersion == kCutOffKernel) {
			for (int i = 0; i < m_num_samples - 1; i++) {
				for (int j = i + 1; j < m_num_samples; j++) {
					if (mv_dis[i][j] < m_dc && mv_dis[i][j] != -1) {
						mv_rho[i] += 1;
						mv_rho[j] += 1;
					}
				}
			}
		}
		else if (LocalDensityVersion == kGaussianKernel) {
			for (int i = 0; i < m_num_samples - 1; i++) {
				for (int j = i + 1; j < m_num_samples; j++) {
					if (mv_dis[i][j] == -1) continue;
					mv_rho[i] = mv_rho[i] + exp(-(mv_dis[i][j] / m_dc) * (mv_dis[i][j] / m_dc));
					mv_rho[j] = mv_rho[j] + exp(-(mv_dis[i][j] / m_dc) * (mv_dis[i][j] / m_dc));
				}
			}
		}
	}

	void DPC::getDistanceToHigherDensity() {
		std::multimap<double, int> sort_rho;
		for (int i = 0; i < m_num_samples; i++)
			sort_rho.insert(std::make_pair(mv_rho[i], i));
		for (auto r_iter = sort_rho.rbegin(); r_iter != sort_rho.rend(); ++r_iter) {
			if (r_iter == sort_rho.rbegin()) {
				mv_delta[r_iter->second] = -1;
				continue;
			}
			for (auto r_iter1 = sort_rho.rbegin(); r_iter1 != r_iter; ++r_iter1) {
				if (mv_dis[r_iter->second][r_iter1->second] < mv_delta[r_iter->second]) {
					mv_delta[r_iter->second] = mv_dis[r_iter->second][r_iter1->second];
					mv_nneigh[r_iter->second] = r_iter1->second;
				}
			}
		}
		double max = -1;
		for (int i = 0; i < m_num_samples; i++) {
			if (max < mv_delta[i]) {
				max = mv_delta[i];
			}
		}
		mv_delta[sort_rho.rbegin()] = max;
	}

	void DPC::getClusterNum() {
		for (int i = 0; i < m_num_samples; i++) {
			if (mv_rho[i] >= m_rho_min && mv_delta[i] >= m_deltamin) {
				mv_cl[i] = m_NCLUST;
				mv_icl.push_back(i);
				++m_NCLUST;
			}
		}
	}

	void DPC::assign() {
		std::multimap<double, int> sort_rho;
		for (int i = 0; i < m_num_samples; i++)
			sort_rho.insert(std::make_pair(mv_rho[i], i));
		for (auto r_iter = sort_rho.rbegin(); r_iter != sort_rho.rend(); ++r_iter) {
			if (mv_cl[r_iter->second] == -1)
				mv_cl[r_iter->second] = mv_cl[mv_nneigh[r_iter->second]];
		}
 	}

	void DPC::getHalo(DPC::Result &result) {
		result.info.resize(m_num_samples);
		for (auto &i : result.info) {
			i.is_center = i.is_core = false;
		}
		result.num_clst = m_NCLUST;

		for (int i = 0; i < m_num_samples; i++)
			mv_halo[i] = mv_cl[i];
		if (m_NCLUST > 1) {
			std::vector<double> bord_rho(m_NCLUST, 0);
			for (int i = 0; i < m_num_samples - 1; i++) {
				for (int j = i + 1; j < m_num_samples; j++) {
					if (mv_cl[i] != mv_cl[j] && mv_dis[i][j] <= m_dc && mv_dis[i][j] != -1) {
						double rho_aver = (mv_rho[i] + mv_rho[j]) / 2.;
						if (rho_aver > bord_rho[mv_cl[i]])
							bord_rho[mv_cl[i]] = rho_aver;
						if (rho_aver > bord_rho[mv_cl[j]])
							bord_rho[mv_cl[j]] = rho_aver;
					}
				}
			}
			for (int i = 0; i < m_num_samples; i++)
				if (mv_rho[i] < bord_rho[mv_cl[i]])
					mv_halo[i] = -1;
		}
		for (int i = 0; i < m_NCLUST; i++) {
			int nc = 0, nh = 0;
			result.info[mv_icl[i]].is_center = true;
			for (int j = 0; j < m_num_samples; j++) {
				if (mv_cl[j] == i) {
					nc++;
					result.info[j].clst_no = i;
				}
				if (mv_halo[j] == i) {
					nh++;
					result.info[j].is_core = true;
				}
			}
			//cout<<endl;
			//cout<<"CLUSTER: "<<i<<" CENTER: "<<mv_icl[i]<<" ELEMENTS: "<<nc<<" CORE: "<<nh<<" HALO: "<<nc-nh<<endl;
		}
	}
	void DPC::autoSetMinRhondMinDelta(double val) {
		//warning: this method need to be investigated
		if (m_num_samples < 2) return;
		std::vector<double> gamma(m_num_samples);
		std::multimap<double, int, std::greater<double>> sortGamma;
		for (int i = 0; i < m_num_samples; ++i) {
			gamma[i] = mv_rho[i] * mv_delta[i];
			sortGamma.insert(std::make_pair(gamma[i], i));
		}
		std::vector<std::pair<int, std::pair<double, double>>> sortPow;
		for (auto i = sortGamma.begin(); i != sortGamma.end(); ++i) {
			sortPow.push_back(std::make_pair(i->second, std::make_pair(log(i->first), log((double)i->second))));
		}
		auto i = sortPow.begin(), j = i + 1;
		double gra1 = (j->second.first - i->second.first) / (j->second.second - i->second.second);
		for (; j != sortPow.end(); ++i, ++j) {
			double gra2 = (j->second.first - i->second.first) / (j->second.second - i->second.second);
			if (fabs(gra2 - gra1) <= val) break;
			gra1 = gra2;
		}
		m_rho_min = (mv_rho[i->first] + mv_rho[j->first]) / 2;
		m_deltamin = (mv_delta[i->first] + mv_delta[j->first]) / 2;
	}

	void DPC::setMinRhoNdMinDelta(double ratio_rho, double ratio_delta) {
		Real minRho, maxRho, minDelta, maxDelta;
		minRho = maxRho = mv_rho[0];
		minDelta = maxDelta = mv_delta[0];
		for (size_t i = 1; i < m_num_samples; i++) {
			if (mv_rho[i] < minRho) minRho = mv_rho[i];
			if (mv_rho[i] > maxRho) maxRho = mv_rho[i];
			if (mv_delta[i] < minDelta) minDelta = mv_delta[i];
			if (mv_delta[i] > maxDelta) maxDelta = mv_delta[i];
		}
		//m_rho_min = ratio * (maxRho - minRho) + minRho;
		//m_deltamin = ratio * (maxDelta - minDelta) + minDelta;
		m_rho_min = ratio_rho * maxRho;
		m_deltamin = ratio_delta * maxDelta;
	}

	void DPC::setMinRhoAndMinDelta(double rho, double delta) {
		m_rho_min = rho;
		m_deltamin = delta;
	}

	void DPC::setDisMatrix(const std::vector<std::vector<double>> &dis) {
		m_num_samples = dis.size();
		mv_rho.resize(m_num_samples);
		for (auto &i : mv_rho) i = 0;
		mv_delta.resize(m_num_samples);
		for (auto &i : mv_delta) i = std::numeric_limits<double>::max();
		mv_cl.resize(m_num_samples);
		for (auto &i : mv_cl) i = -1;
		mv_nneigh.resize(m_num_samples);
		for (auto &i : mv_nneigh) i = -1;
		mv_halo.resize(m_num_samples);
		for (auto &i : mv_halo) i = 0;
		mv_dis = dis;
		m_NCLUST = 0;
		m_dc = 0.0;
		m_rho_min = 1;
		m_deltamin = 0.1;
	}

	void DPC::clustering(DPC::Result &result) {
		if (m_num_samples <= 0) return;
		getdc();
		getLocalDensity();
		getDistanceToHigherDensity();
		//set rhomin and deltamin here
		autoSetMinRhondMinDelta(0.1);
		getClusterNum();
		assign();
		getHalo(result);
	}
	
	void DPC::clustering() {
		if (m_num_samples <= 0) return;
		getdc();
		getLocalDensity();
		getDistanceToHigherDensity();
		//set rhomin and deltamin here
		autoSetMinRhondMinDelta(0.1);
		getClusterNum();
		assign();
	}

	void DPC::clustering(double ratio_rho, double ratio_delta) {
		if (m_num_samples <= 0) return;
		getdc();
		getLocalDensity(kCutOffKernel);
		getDistanceToHigherDensity();
		setMinRhoNdMinDelta(ratio_rho, ratio_delta);
		getClusterNum();
		assign();
	}

	std::vector<std::vector<size_t>> DPC::clusters() const {
		std::vector<std::vector<size_t>> clusters(m_NCLUST);
		for (size_t i = 0; i < m_num_samples; ++i) {
			clusters[mv_cl[i]].push_back(i);
		}
		return clusters;
	}
}