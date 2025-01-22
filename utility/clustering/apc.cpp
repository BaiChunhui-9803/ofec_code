#include "apc.h"

namespace ofec {
	void APC::clustering(Problem *pro) {
		updateSimilarity(pro);
		size_t num_iter_stagnant = 0;
		for (size_t iter = 0; iter < m_Mits; iter++) {
			updateResponsibility();
			updateAvailability();
			std::vector<bool> last_is_exemplar(m_is_exemplar);
			size_t num_exemplars = 0;
			for (size_t i = 0; i < m_N; ++i) {
				m_is_exemplar[i] = m_responsibility[i][i] + m_availability[i][i] > 0;
				if (m_is_exemplar[i])
					num_exemplars++;
			}
			if (num_exemplars == 0) {
				num_iter_stagnant = 0;
				continue;
			}
			bool same = true;
			for (size_t i = 0; i < m_N; ++i) {
				if (m_is_exemplar[i] != last_is_exemplar[i]) {
					same = false;
					break;
				}
			}
			if (same)
				num_iter_stagnant++;
			else
				num_iter_stagnant = 0;
			if (num_iter_stagnant == m_Cits)
				break;
		}
		updateClusters();
	}

	void APC::updateSimilarity(Problem *pro) {
		std::vector<Real> all_similarities;
		for (size_t i = 0; i < m_N - 1; i++) {
			for (size_t j = i + 1; j < m_N; j++) {
				m_similarity[i][j] = m_similarity[j][i] = -pow(pro->variableDistance(
					m_data[j]->variableBase(), m_data[i]->variableBase()), 2);
				all_similarities.emplace_back(m_similarity[i][j]);
			}
		}
		std::sort(all_similarities.begin(), all_similarities.end());
		Real median;
		size_t size = all_similarities.size();
		if (size % 2 == 0)
			median = (all_similarities[size / 2] + all_similarities[size / 2 - 1]) / 2;
		else
			median = all_similarities[size / 2];
		for (size_t i = 0; i < m_N; i++) {
			m_similarity[i][i] = median;
		}
	}

	void APC::updateResponsibility() {
		/* original version */
		//for (size_t i = 0; i < m_N; i++) {
		//	for (size_t k = 0; k < m_N; k++) {
		//		Real max = -std::numeric_limits<Real>::max();
		//		for (size_t kk = 0; kk < m_N; ++kk) {
		//			if (kk == k)
		//				continue;
		//			if (max < m_availability[i][kk] + m_similarity[i][kk])
		//				max = m_availability[i][kk] + m_similarity[i][kk];
		//		}
		//		m_responsibility[i][k] = m_lambda * m_responsibility[i][k] + (1 - m_lambda) * (m_similarity[i][k] - max);
		//	}
		//}

		/* optimized version */
		std::vector<Real> max(m_N), sec(m_N);
		for (size_t i = 0; i < m_N; i++) {
			max[i] = std::max(m_availability[i][0] + m_similarity[i][0], m_availability[i][1] + m_similarity[i][1]);
			sec[i] = std::min(m_availability[i][0] + m_similarity[i][0], m_availability[i][1] + m_similarity[i][1]);
			for (size_t k = 2; k < m_N; k++) {
				Real val = m_availability[i][k] + m_similarity[i][k];
				if (val < sec[i])
					continue;
				else if (val < max[i])
					sec[i] = val;
				else {
					sec[i] = max[i]; max[i] = val;
				}
			}
		}
		for (size_t i = 0; i < m_N; i++) {
			for (size_t k = 0; k < m_N; k++) {
				if (m_availability[i][k] + m_similarity[i][k] == max[i])
					m_responsibility[i][k] = m_lambda * m_responsibility[i][k] + (1 - m_lambda) * (m_similarity[i][k] - sec[i]);
				else
					m_responsibility[i][k] = m_lambda * m_responsibility[i][k] + (1 - m_lambda) * (m_similarity[i][k] - max[i]);
			}
		}
	}

	void APC::updateAvailability() {
		/* original version */
		//for (size_t i = 0; i < m_N; ++i) {
		//	for (size_t k = 0; k < m_N; ++k) {
		//		if (i == k) {
		//			Real sum = 0;
		//			for (size_t ii = 0; ii < m_N; ++ii) {
		//				if (ii == k)
		//					continue;
		//				sum += std::max<Real>(0, m_responsibility[ii][k]);
		//			}
		//			m_availability[i][k] = m_lambda * m_availability[i][k] + (1 - m_lambda) * sum;
		//		}
		//		else {
		//			Real sum = 0;
		//			for (size_t ii = 0; ii < m_N; ++ii) {
		//				if (ii == i || ii == k)
		//					continue;
		//				sum += std::max<Real>(0, m_responsibility[ii][k]);
		//			}
		//			m_availability[i][k] = m_lambda * m_availability[i][k] + (1 - m_lambda) * std::min<Real>(0, m_responsibility[k][k] + sum);
		//		}
		//	}
		//}

		/* optimized version  */
		std::vector<Real> sum(m_N, 0);
		for (size_t k = 0; k < m_N; k++) {
			for (size_t ii = 0; ii < m_N; ii++) {
				sum[k] += std::max<Real>(0, m_responsibility[ii][k]);
			}
		}
		for (size_t i = 0; i < m_N; i++) {
			for (size_t k = 0; k < m_N; k++) {
				if (i == k)
					m_availability[i][k] = m_lambda * m_availability[i][k] + (1 - m_lambda) * (sum[k] - std::max<Real>(0, m_responsibility[k][k]));
				else
					m_availability[i][k] = m_lambda * m_availability[i][k] + (1 - m_lambda) * (std::min<Real>(0, m_responsibility[k][k] + (sum[k] - std::max<Real>(0, m_responsibility[i][k]) - std::max<Real>(0, m_responsibility[k][k]))));
			}
		}
	}

	void APC::updateClusters() {
		std::vector<size_t> vec_exemplars;
		for (size_t i = 0; i < m_N; i++) {
			if (m_is_exemplar[i]) {
				vec_exemplars.emplace_back(i);
				m_clusters.emplace_back(std::vector<size_t>({ i }));
			}
		}
		if (m_clusters.empty()) {
			std::vector<size_t> all_inds(m_data.size());
			std::iota(all_inds.begin(), all_inds.end(), 0);
			m_clusters.emplace_back(all_inds);
			return;
			//			std::ofstream out_file("C:/Users/WJC/Desktop/data.csv");
			//			for (size_t i = 0; i < m_N; i++) {
			//				for (Real val : m_data[i]->variable().vect())
			//					out_file << val << ",";
			//				out_file << std::endl;
			//			}
						//throw("error at updateClusters");
		}
		for (size_t i = 0; i < m_N; i++) {
			if (m_is_exemplar[i])
				continue;
			size_t idx = 0;
			for (size_t k = 1; k < vec_exemplars.size(); ++k) {
				if (m_similarity[i][vec_exemplars[k]] > m_similarity[i][vec_exemplars[idx]])
					idx = k;
			}
			m_clusters[idx].emplace_back(i);
		}
	}
}