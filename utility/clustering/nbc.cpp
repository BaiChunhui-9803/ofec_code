#include "nbc.h"
#include <map>
#include "../functional.h"
#include "../kd-tree/data_adaptor.h"
#include "../../core/problem/continuous/continuous.h"
#include "../../core/problem/solution.h"
#include "../functional.h"
#include "../../core/environment/environment.h"

namespace ofec {
	NBC::NBC(Real phi, UpdateNBD update_nbd, CalThrshld cal_thrdshld, bool use_rule2): 
		m_phi(phi), m_N(0), m_update_nbd(update_nbd), m_cal_thrshld(cal_thrdshld),
		m_mean(0), m_stddev(0), m_use_rule2(use_rule2) {}

	void NBC::setData(const std::vector<const SolutionBase*> &sols) {
		m_N = sols.size();
		m_graph.resize(m_N);
		m_data = sols;
		if (!m_clusters.empty())
			m_clusters.clear();
	}

	void NBC::clustering(Environment *env) {
		switch (m_update_nbd) {
		case kByDistMat:
			updateNbDistByDistMat(env);
			break;
		case kByKDTree:
			if (CAST_CONOP(env->problem()) == nullptr) {
				throw Exception("NBD cannot be calculated by k-d tree search when the problem is not continuous");
			}			
			updateNbDistByKDTree(env);
			break;
		}
		
		// RULE1
		cutEdgesInGraph();
		// RULE2
		if (m_use_rule2 && CAST_CONOP(env->problem()) != nullptr) {
			std::vector<std::list<size_t>> incoming(m_N);
			for (size_t i = 0; i < m_N; ++i) {
				if (m_graph[i] != -1) {
					incoming[m_graph[i]].push_back(i);
				}
			}
			// b(S, D) = (-4.69 * 10-4 * D2 + 0.0263 * D + 3.66/D - 0.457) * log10(S) + 7.51e - 4 * D2 - 0.0421 * D - 2.26 / D + 1.83
			Real D = env->problem()->numberVariables();
			Real D2 = pow(D, 2);
			Real S = m_N;
			Real b = (-4.69 * 10 - 4 * D2 + 0.0263 * D + 3.66 / D - 0.457) * log10(S) + 7.51e-4 * D2 - 0.0421 * D - 2.26 / D + 1.83;
			for (size_t i = 0; i < m_N; ++i) {
				if (m_graph[i] != -1 && incoming[i].size() >= 3) {
					std::vector<Real> nbds;
					for (size_t k : incoming[i]) {
						nbds.push_back(m_nb_dis_2[k]);
					}
					std::nth_element(nbds.begin(), nbds.begin() + nbds.size() / 2, nbds.end());
					Real median_nbd = *(nbds.begin() + nbds.size() / 2);
					if (m_nb_dis_2[i] / median_nbd > b) {
						m_graph[i] = -1;
					}
				}
			}
		}

		updateClusters();
	}

    void NBC::clustering(size_t num, Environment *env) {
        if(num < 2) {
            m_graph[0] = -1;
            m_graph[1] = 0;
        }
        else {
            switch (m_update_nbd) {
                case kByDistMat:
                    updateNbDistByDistMat(env);
                    break;
                case kByKDTree:
                    updateNbDistByKDTree(env);
                    break;
            }
			cutEdgesInGraph(num);
        }
        updateClusters();
    }

	void NBC::updateNbDistByDistMat(Environment *env) {
		/* update the distance matrix */
		if (m_distance.size() != m_N) {
			m_distance.resize(m_N);
			for (auto &row : m_distance)
				row.assign(m_N, 0);
		}
		for (size_t i = 0; i < m_N; i++) {
			for (size_t k = i + 1; k < m_N; k++) {
				m_distance[i][k] = m_distance[k][i] = m_data[i]->variableDistance(*m_data[k], env);
			}
		}

		/* update the nearest-better distance */
		m_graph.assign(m_N, -1);
		for (size_t i = 0; i < m_N; i++) {
			std::map<Real, size_t> seq_dis;           // the default hehavior is in ascending order
			for (size_t k = 0; k < m_N; k++) {
				if (i != k)
					seq_dis.emplace(std::make_pair(m_distance[i][k], k));
			}
			for (const auto &p : seq_dis) {
				if (dominate(*m_data[p.second],*m_data[i], env->problem()->optimizeMode())
					|| (env->problem()->same(m_data[p.second]->variableBase(), 
						m_data[i]->variableBase()) && p.second < i))
				{
					m_graph[i] = p.second;
					break;
				}
			}
		}

		m_nb_dis_1.assign(m_N - 1, -1);
		m_nb_dis_2.assign(m_N, -1);
		for (size_t i = 0, k = 0; i < m_N; ++i) {
			if (m_graph[i] != -1)
				m_nb_dis_2[i] = m_nb_dis_1[k++] = m_distance[i][m_graph[i]];
		}
	}

	void NBC::updateNbDistByKDTree(Environment *env) {
		if (CAST_CONOP(env->problem()) == nullptr)
			throw Exception("Calculation of NBDs by k-D tree is only feasible in continuous search space");

		std::vector<std::vector<Real>> samples(m_data.size());
		for (size_t i = 0; i < m_data.size(); ++i)
			samples[i] = dynamic_cast<const Solution<>*>(m_data[i])->variable().vector();

		typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real>  MyKdTree;
		size_t num_vars = CAST_CONOP(env->problem())->numberVariables();
		MyKdTree mat_index(num_vars, samples);
		size_t num_nn = 3;
		std::set<size_t> fnd;
		m_graph.assign(m_N, -1);
		m_nb_dis_1.assign(m_N - 1, -1);
		m_nb_dis_2.assign(m_N, -1);
		while (fnd.size() < samples.size() - 1) {
			std::vector<size_t> ret_indexes(num_nn);
			std::vector<Real> out_dists_sqr(num_nn);
			nanoflann::KNNResultSet<Real> result_set(num_nn);
			for (size_t i = 0; i < samples.size(); ++i) {
				if (fnd.count(i))
					continue;
				result_set.init(&ret_indexes[0], &out_dists_sqr[0]);
				mat_index.index->findNeighbors(result_set, &samples[i][0], nanoflann::SearchParams(10));
				for (size_t n = 0; n < num_nn; ++n) {
					if (dominate(*m_data[ret_indexes[n]], *m_data[i], env->problem()->optimizeMode())
						|| (m_data[ret_indexes[n]]->objectiveDistance(*m_data[i]) == 0 && ret_indexes[n] < i))
					{
                        m_nb_dis_2[i] = m_nb_dis_1[fnd.size()] = sqrt(out_dists_sqr[n]);
                        m_graph[i] = ret_indexes[n];
                        fnd.insert(i);
                        break;
                    }
				}
			}
			num_nn = 1 + (num_nn - 1) * 2;
			if (num_nn > samples.size())
				num_nn = samples.size();
		}
	}

	void NBC::cutEdgesInGraph() {
		Real thrshld;

		switch (m_cal_thrshld) {
			case kByOutlier: {
				if (m_nb_dis_1.size() < 100) {
					auto nb_dis_1 = m_nb_dis_1;
					std::vector<Real> nb_dis_2;
					while (!nb_dis_1.empty()) {
						Real min = nb_dis_1[0];
						size_t id_min = 0;
						for (size_t i = 0; i < nb_dis_1.size(); ++i) {
							if (min > nb_dis_1[i]) {
								id_min = i;
								min = nb_dis_1[i];
							}
						}
						nb_dis_2.push_back(min);
						calMeanAndStd(nb_dis_2, m_mean, m_stddev);
						thrshld = m_mean + m_phi * m_stddev;
						if (min > thrshld)
							break;
						else {
							nb_dis_1[id_min] = nb_dis_1.back();
							nb_dis_1.pop_back();
						}
					}
				}
				else {
					auto nb_dis = m_nb_dis_1;
					calMeanAndStd(nb_dis, m_mean, m_stddev);
					while (true) {
						thrshld = m_mean + m_phi * m_stddev;
						Real max = nb_dis[0];
						size_t id_max = 0;
						for (size_t i = 0; i < nb_dis.size(); ++i) {
							if (max < nb_dis[i]) {
								id_max = i;
								max = nb_dis[i];
							}
						}
						if (max <= thrshld)
							break;
						else {
							nb_dis[id_max] = nb_dis.back();
							nb_dis.pop_back();
							calMeanAndStd(nb_dis, m_mean, m_stddev);
						}
					}
				}
				break;
			}
			case kByMean: default: {
				calMeanAndStd(m_nb_dis_1, m_mean, m_stddev);
				thrshld = m_phi * m_mean;
			}
		}

		for (size_t i = 0; i < m_N; i++) {
			if (m_graph[i] != -1 && m_nb_dis_2[i] > thrshld) {
				m_graph[i] = -1;
			}
		}
	}

    void NBC::cutEdgesInGraph(size_t num) {
        // sort
        auto dis = m_nb_dis_1;
        std::nth_element(dis.begin(), dis.end() - num + 1, dis.end());
        calMeanAndStd(m_nb_dis_1, m_mean, m_stddev);

        for (size_t i = 0; i < m_N; i++) {            // cut edges too long
			if (m_graph[i] != -1 && m_nb_dis_2[i] >= *(dis.end() - num + 1)) {
				m_graph[i] = -1;
			}
        }
    }

	void NBC::cutEdgesInGraph(Real thrshld) {
		for (size_t i = 0; i < m_N; i++) {            // cut edges too long
			if (m_graph[i] != -1 && m_nb_dis_2[i] > thrshld) {
				m_graph[i] = -1;
			}
		}
	}

	void NBC::cutEdgesInGraph(std::vector<size_t> &ids_centers) {
		for (size_t &id_center : ids_centers) {
			while (m_nb_dis_2[id_center] != -1 && m_nb_dis_2[id_center] == 0) {
				id_center = m_graph[id_center];
			}
			m_graph[id_center] = -1;
		}
	}

	void NBC::updateClusters() {
		m_cluster_centers.clear();
		for (size_t i = 0; i < m_N; i++) {
			if (m_graph[i] == -1) {
				m_cluster_centers.push_back(i);
			}
		}

		std::vector<int> idxs_clu(m_N, -1);
		size_t  idx_clu(0), idx_cur_ind;
		for (size_t i = 0; i < m_N; i++) {
			if (idxs_clu[i] == -1) {
				std::list<size_t> path;
				idx_cur_ind = i;
				path.emplace_back(idx_cur_ind);
				while (m_graph[idx_cur_ind] != -1) {
					idx_cur_ind = m_graph[idx_cur_ind];
					path.emplace_back(idx_cur_ind);
				}
				if (idxs_clu[idx_cur_ind] == -1)
					idxs_clu[idx_cur_ind] = idx_clu++;
				for (size_t idx_ind : path) {
					idxs_clu[idx_ind] = idxs_clu[idx_cur_ind];
					m_graph[idx_ind] = -1;
				}
			}
		}
		m_clusters.clear();
		m_clusters.resize(idx_clu);
		for (size_t i = 0; i < m_N; ++i) {
			m_clusters[idxs_clu[i]].emplace_back(i);
		}
	}

}