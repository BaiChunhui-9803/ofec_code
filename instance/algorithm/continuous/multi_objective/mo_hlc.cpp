#include "mo_hlc.h"
#include "../../../../core/problem/continuous/continuous.h"

namespace ofec {
	MO_HLC::MO_HLC(size_t dim_var, size_t num_obj):m_split_tree_MOP(std::make_unique<KDTree_MOP>()),\
		m_selection_tree_MOP(std::make_unique<KDTree_MOP>()), m_obj_tree_MOP(std::make_unique<KDTree_MOP>()) {
		//std::cout << "call this constructor" << std::endl;
	}

	void MO_HLC::initialVarSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num) {
		/*initialize subspace of MOP*/
		for (int i = 0; i < subspace_num; ++i) {
			m_subspace_info.emplace_back(new SubspaceInfo_MOP);
			//m_basin_info(new BasinInfo_MOP);//在子空间聚类之后更新
		}
		
		m_split_tree_MOP->setInitBox(boundary);
		m_split_tree_MOP->inputRatioData(std::vector<Real>(subspace_num, 1.0 / subspace_num));
		m_split_tree_MOP->buildIndex();
		for (size_t i = 0; i < m_subspace_info.size(); ++i) {
			m_split_tree_MOP->findNeighbor(i, m_subspace_info[i]->m_sub_neighbors);
			m_subspace_info[i]->subspace_inx = i;
		}

		Real volume = 1.;
		for (size_t i = 0; i < boundary.size(); ++i) {
			volume *= (boundary[i].second - boundary[i].first);
		}
		setSpaceVoume(volume);
	}

	void MO_HLC::initialVarSelectionSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num) {
		/*initialize variable space selection subspaces of MOP*/
		m_selection_tree_MOP->setInitBox(boundary);
		m_selection_tree_MOP->inputRatioData(std::vector<Real>(subspace_num, 1.0 / subspace_num));
		m_selection_tree_MOP->buildIndex();
	}

	void MO_HLC::initialObjSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num) {
		/*initialize subspace of MOP*/
		m_obj_space_info.clear();
		for (int i = 0; i < subspace_num; ++i) {
			m_obj_space_info.emplace_back(new ObjRegionInfo);
		}

		m_obj_tree_MOP->setInitBox(boundary);
		m_obj_tree_MOP->inputRatioData(std::vector<Real>(subspace_num, 1.0 / subspace_num));
		m_obj_tree_MOP->buildIndex();
		for (size_t i = 0; i < m_obj_space_info.size(); ++i) {
			m_obj_tree_MOP->findNeighbor(i, m_obj_space_info[i]->m_obj_neighbors);
			m_obj_space_info[i]->subspace_inx = i;
		}
	}

	/*采样解*/
	void MO_HLC::getSampling(const Solution<>& sample_point) {
		/* Update statistics */
		size_t idx_region = m_split_tree_MOP->getRegionIdx(sample_point.variable().vect());
		//m_subspace_info[idx_region].m_current_size++;
		//m_subspace_info[idx_region].m_total_size++;
		//依次与子空间内各层个体进行支配关系比较
		//if (!m_subspace_info[idx_region].best_sol || sample_point.dominate(*m_region_info[idx_region].best_sol))//采用多目标的比较方式
		//	m_region_info[idx_region].stat[0] = sample_point.objective(0);
		////m_region_info[idx_region].stat[0] = (m_region_info[idx_region].stat[0] * m_region_info[idx_region].stat[1] + sample_point.objective(0)) / (m_region_info[idx_region].stat[1] + 1);
		//m_region_info[idx_region].stat[1]++;
		//m_his_inds.emplace_back(sample_point);
		//m_region_info[idx_region].sols.emplace_back(&m_his_inds.back());
		//if (!m_region_info[idx_region].best_sol || m_his_inds.back().dominate(*m_region_info[idx_region].best_sol)) {
		//	m_region_info[idx_region].best_sol = &m_his_inds.back();
		//}
		//m_num_sample++;

		//if (m_worst_ind == nullptr || m_worst_ind->dominate(m_his_inds.back()))
		//	m_worst_ind = &m_his_inds.back();
	}

	void MO_HLC::initialBasin(size_t num_basin) {
		//size_t num = getClusters();
		for (size_t i = 0; i < num_basin; ++i) {
			m_basin_info.emplace_back(new BasinInfo_MOP);
		}
	}

	//每次聚完一个类
	//void MO_HLC::rank_clustering() {
	//	m_clusters.clear();
	//	std::vector<size_t> cluster;
	//	std::vector<size_t> temp_idx(m_num_regions, 0);
	//	//先找出排序值最靠前的子空间
	//	//std::vector<size_t> front_subspace;
	//	while (std::find(temp_idx.begin(), temp_idx.end(), 0) != temp_idx.end()) {
	//		//先得到现有子空间最好的排序值的索引
	//		Real rank = INT16_MAX;//存放最靠前的rank值
	//		size_t idx = 0;
	//		for (size_t i = 0; i < m_subspace_info.size(); ++i) {
	//			if (temp_idx[i] == 0 && m_subspace_info[i].m_best_rank < rank) {
	//				rank = m_subspace_info[i].m_best_rank;
	//				idx = i;
	//			}
	//		}
	//		std::vector<size_t> current_idx;
	//		current_idx.push_back(idx);
	//		cluster.push_back(idx);
	//		temp_idx[idx] = 1;
	//		while (!current_idx.empty()) {
	//			std::vector<size_t> temp_current;
	//			for (size_t j = 0; j < current_idx.size(); ++j) {
	//				size_t temp_rank = m_subspace_info[current_idx[j]].m_best_rank;
	//				auto temp_neighbor = m_neighbor[current_idx[j]];
	//				std::vector<size_t> sub;
	//				for (auto& k : temp_neighbor) {
	//					if (temp_idx[k] == 0) {
	//						sub.push_back(k);
	//					}
	//				}
	//				for (size_t p = 0; p < sub.size(); ++p) {
	//					if (m_subspace_info[sub[p]].m_best_rank >= temp_rank&& m_subspace_info[sub[p]].m_best_rank!=0) {
	//						cluster.push_back(sub[p]);
	//						temp_idx[sub[p]] = 1;
	//						temp_current.push_back(sub[p]);
	//					}
	//				}
	//			}
	//			//更新current_idx
	//			if (!temp_current.empty())
	//				current_idx = temp_current;
	//			else
	//				break;
	//		}
	//		m_clusters.emplace_back(cluster);
	//		cluster.clear();
	//	}
	//}

	//每次聚类只加入一层，允许类间有重叠
	void MO_HLC::rankClustering() {
		m_clusters_MOP.clear();
		//std::vector<size_t> cluster;
		std::vector<size_t> flag_idx(m_subspace_info.size(), 0);
		//先找出排序值最靠前的子空间
		while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
			//先得到现有子空间最好的排序值的索引
			int rank = INT16_MAX;//存放最靠前的rank值
			std::vector<size_t> best_idx;//具有最好rank值的子空间索引
			for (size_t i = 0; i < m_subspace_info.size(); ++i) {
				if (flag_idx[i] == 0 && m_subspace_info[i]->m_best_rank < rank) {
					rank = m_subspace_info[i]->m_best_rank;
				}
			}
			for (size_t i = 0; i < m_subspace_info.size(); ++i) {
				if (flag_idx[i] == 0 && m_subspace_info[i]->m_best_rank == rank)
					best_idx.push_back(i);
			}
			//对这些最靠前的子空间进行聚类
			std::vector<std::vector<size_t>> cluster_idx;
			cluster_idx = clusterFrontSpace(best_idx, false);

			for (size_t i = 0; i < cluster_idx.size(); ++i) {
				for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
					flag_idx[cluster_idx[i][j]] = 1;
				}
			}

			bool gradual_cluster = true;//是否逐层扩散
			bool strict_cluster = true;//对于每个子空间，只有其未聚类的邻域排序值全部小于该子空间，才加入到聚类
			bool permit_overlap = false;
			if (gradual_cluster) {
				auto current = cluster_idx;
				while (!current.empty()) {
					std::vector<std::vector<size_t>> temp_current;//所有类的当前比较的子空间
					for (size_t i = 0; i < current.size(); ++i) {
						std::vector<size_t> temp1;//每一个类的第一层邻域
						if (current[i].size() > 1 || current[i][0] != INT16_MAX) {
							for (size_t j = 0; j < current[i].size(); ++j) {
								size_t temp_rank = m_subspace_info[current[i][j]]->m_best_rank;
								auto temp_neighbor = m_subspace_info[current[i][j]]->m_sub_neighbors;
								std::vector<size_t> sub;
								for (auto& k : temp_neighbor) {
									if (permit_overlap) {//是否允许聚类的子空间重叠
										sub.push_back(k);
									}
									else {
										if (flag_idx[k] == 0) {
											sub.push_back(k);
										}
									}
								}
								if (strict_cluster) {
									std::vector<Real> temp_space;
									bool add = true;
									for (size_t p = 0; p < sub.size(); ++p) {
										if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
											if (m_subspace_info[sub[p]]->m_best_rank > temp_rank ) {
												temp_space.push_back(sub[p]);
											}
											else {
												add = false;
												break;
											}
										}
									}
									if (add) {
										for (auto& k : temp_space) {
											cluster_idx[i].push_back(k);
											temp1.push_back(k);
											flag_idx[k] = 1;
										}
									}
								}
								else {
									for (size_t p = 0; p < sub.size(); ++p) {
										if (m_subspace_info[sub[p]]->m_best_rank > temp_rank ) {
											if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
												cluster_idx[i].push_back(sub[p]);
												temp1.push_back(sub[p]);
												flag_idx[sub[p]] = 1;
											}
										}
									}
								}
							}
						}
						if (temp1.empty()) {
							temp1.push_back(INT16_MAX);
						}
						temp_current.push_back(temp1);
					}
					bool over = false;
					for (size_t i = 0; i < temp_current.size(); ++i) {
						if (temp_current[i][0] != INT16_MAX) {
							break;
						}
						else if (i == temp_current.size() - 1) {
							over = true;
						}
					}
					if (over) {
						current.clear();
					}
					else {
						current = temp_current;
						temp_current.clear();
					}
				}
			}
			else {
				for (size_t i = 0; i < cluster_idx.size(); ++i) {
					std::vector<size_t> cur_compare_space = cluster_idx[i];
					while (!empty(cur_compare_space)) {
						auto compare_space = cur_compare_space;
						cur_compare_space.clear();
						for (size_t j = 0; j < compare_space.size(); ++j) {
							auto neighbors = m_subspace_info[compare_space[j]]->m_sub_neighbors;
							std::vector<size_t> sub;
							for (auto& k : neighbors) {
								if (permit_overlap) {//是否允许聚类的子空间重叠
									sub.push_back(k);
								}
								else {
									if (flag_idx[k] == 0) {
										sub.push_back(k);
									}
								}
							}
							if (strict_cluster) {
								std::vector<Real> temp_space;
								bool add = true;
								for (size_t p = 0; p < sub.size(); ++p) {
									if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
										if (m_subspace_info[sub[p]]->m_best_rank > m_subspace_info[compare_space[j]]->m_best_rank) {
											temp_space.push_back(sub[p]);
											cluster_idx[i].push_back(sub[p]);
											cur_compare_space.push_back(sub[p]);
											flag_idx[sub[p]] = 1;
										}
										else {
											add = false;
											break;
										}
									}
								}
								if (add) {
									for (auto& k : temp_space) {
										cluster_idx[i].push_back(k);
										cur_compare_space.push_back(k);
										flag_idx[k] = 1;
									}
								}
							}
							else {
								for (size_t p = 0; p < sub.size(); ++p) {
									if (m_subspace_info[sub[p]]->m_best_rank > m_subspace_info[compare_space[j]]->m_best_rank) {
										if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
											cluster_idx[i].push_back(sub[p]);
											cur_compare_space.push_back(sub[p]);
											flag_idx[sub[p]] = 1;
										}
									}
								}
							}
						}
					}
				}
			}
			std::vector<std::vector<size_t>> temp_clu;
			for (size_t i = 0; i < cluster_idx.size(); ++i) {
				temp_clu.emplace_back(cluster_idx[i]);
			}
			m_clusters_MOP.emplace_back(temp_clu);
		}
	}

	//根据潜力值聚类
	std::vector<std::vector<std::vector<size_t>>> MO_HLC::clusterSubspace(const std::vector<size_t>& subspaces, size_t total_num_spaces,bool neighbor_flag) {//根据潜力值聚类
		m_clusters_MOP.clear();
		std::vector<std::vector<std::vector<size_t>>> cluster;
		//std::vector<size_t> flag_idx(total_num_spaces, 0);
		//for (size_t i = 0; i < total_num_spaces; ++i) {
		//	if (std::find(subspaces.begin(), subspaces.end(), i) == subspaces.end()) {
		//		flag_idx[i] = 1;
		//	}
		//}
		//while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
		//	//先得到现有子空间最好指标值的索引
		//	Real potential = -1.*INT16_MAX;//存放最好的指标值
		//	std::vector<size_t> best_idx;//具有最好指标值的子空间索引
		//	for (size_t i = 0; i < total_num_spaces; ++i) {
		//		if (flag_idx[i] == 0 && m_subspace_info[i]->m_potential > potential) {
		//			potential = m_subspace_info[i]->m_potential;
		//		}
		//	}
		//	for (size_t i = 0; i < total_num_spaces; ++i) {
		//		if (flag_idx[i] == 0 && m_subspace_info[i]->m_potential == potential)
		//			best_idx.push_back(i);
		//	}
		//	//对这些最靠前的子空间进行聚类
		//	std::vector<std::vector<size_t>> cluster_idx;
		//	std::vector<size_t> flag(best_idx.size(), 0);
		//	cluster_idx = clusterFrontSpace(best_idx,neighbor_flag);
		//	
		//	//标记已聚类子空间
		//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
		//		for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
		//			flag_idx[cluster_idx[i][j]] = 1;
		//		}
		//	}

		//	//然后根据邻域潜力值向外依次逐层扩散
		//	auto current = cluster_idx;
		//	while (!current.empty()) {
		//		std::vector<std::vector<size_t>> temp_current;
		//		for (size_t i = 0; i < current.size(); ++i) {
		//			std::vector<size_t> temp1;
		//			for (size_t j = 0; j < current[i].size(); ++j) {
		//				if (current[i][j]!=INT16_MAX) {
		//					Real temp_potential = m_subspace_info[current[i][j]]->m_potential;
		//					auto temp_neighbor = m_subspace_info[current[i][j]]->m_sub_neighbors;
		//					std::vector<size_t> sub;
		//					for (auto& k : temp_neighbor) {
		//						/*if (flag_idx[k] == 0) {
		//							sub.push_back(k);
		//						}*/
		//						sub.push_back(k);
		//					}
		//					for (size_t p = 0; p < sub.size(); ++p) {
		//						if (std::fabs(m_subspace_info[sub[p]]->m_potential - temp_potential)/ temp_potential <=0.1) {//模糊聚类
		//							cluster_idx[i].push_back(sub[p]);
		//							temp1.push_back(sub[p]);
		//							flag_idx[sub[p]] = 1;
		//						}
		//					}
		//				}
		//			}
		//			if (temp1.empty()) {
		//				temp1.push_back(INT16_MAX);
		//			}
		//			temp_current.push_back(temp1);
		//		}
		//		bool over = false;
		//		for (size_t i = 0; i < temp_current.size(); ++i) {
		//			if (temp_current[i][0] !=INT16_MAX) {
		//				break;
		//			}
		//			else if (i == temp_current.size() - 1) {
		//				over = true;
		//			}
		//		}
		//		if (over) {
		//			current.clear();
		//		}
		//		else {
		//			current = temp_current;
		//		}
		//	}
		//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
		//		cluster.push_back(cluster_idx[i]);
		//	}
		//}
		m_clusters_MOP = cluster;
		return cluster;
	}

	std::vector<std::vector<size_t>> MO_HLC::clusterFrontSpace(const std::vector<size_t>& frontspace,bool neighbor_flag) {
		//进行最优子空间的聚类
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space = i;
					head_cluster.push_back(frontspace[begin_space]);
					select_flag[i] = 1;
					break;
				}
			}
			for (size_t i = 0; i < select_flag.size(); ++i) {
				count += select_flag[i];
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
			auto temp_cluster = head_cluster;
			while (count< select_flag.size()) {
				std::vector<size_t> temp;
				for (size_t j = 0; j < temp_cluster.size(); ++j) {
					size_t inx = temp_cluster[j];
					std::list<size_t> neighbors;
					subspaceTree().findNeighbor(inx, neighbors);
					for (size_t k = 0; k < frontspace.size(); ++k) {
						if (select_flag[k] == 0) {
							if (std::find(neighbors.begin(), neighbors.end(), frontspace[k]) != neighbors.end()) {
								head_cluster.push_back(frontspace[k]);
								temp.push_back(frontspace[k]);
								select_flag[k] = 1;
								count++;
							}
						}
					}
				}
				if (temp.empty()) {
					clustered.emplace_back(head_cluster);
					break;
				}
				else {
					temp_cluster = temp;
				}
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
		}
		if (neighbor_flag) {
			//再加入各类的一级邻域
			for (size_t i = 0; i < clustered.size(); ++i) {
				std::vector<size_t> m_neighbors;
				for (size_t j = 0; j < clustered[i].size(); ++j) {
					size_t inx = clustered[i][j];
					std::list<size_t> neighbors;
					subspaceTree().findNeighbor(inx, neighbors);
					for (auto k : neighbors) {
						if (m_neighbors.empty() || std::find(m_neighbors.begin(), m_neighbors.end(), k) == m_neighbors.end()) {
							m_neighbors.push_back(k);
						}
					}
				}
				for (auto j : m_neighbors) {
					if (std::find(clustered[i].begin(), clustered[i].end(), j) == clustered[i].end()) {
						clustered[i].push_back(j);
					}
				}
			}
		}
		m_front_clusters = clustered;
		return clustered;
	}

	void MO_HLC::findClusterCenterSsp() {

	}

	void MO_HLC::updateSubspaceInfo(Population<Solution<>> &pop,Problem *pro,Random *rnd) {
		Real total_volume = getSpaceVolume();
		Real min_volume = total_volume;
		//更新子空间邻域关系
		for (size_t i = 0; i < numSubspace(); ++i) {
			subspaceTree().findNeighbor(i, getSubspaceInfo(i).m_sub_neighbors);
			auto temp_volume = subspaceTree().getBoxVolume(i);
			if (temp_volume < min_volume) {
				min_volume = temp_volume;
			}
		}
		std::map<size_t, std::vector<size_t>> pop2space;
		for (size_t i = 0; i < pop.size(); ++i) {
			auto var = pop[i].variable().vect();
			auto idx = subspaceTree().getRegionIdx(var);
			if (pop2space[idx].empty()) {
				std::vector<size_t> temp_inx;
				pop2space.insert(std::make_pair(idx, temp_inx));
			}
			pop2space[idx].emplace_back(i);
		}
		//更新子空间前沿解和代表解
		for (auto& sub : pop2space) {
			size_t idx = sub.first;
			Population<Solution<>> temp_pop;
			for (size_t i = 0; i < sub.second.size(); ++i) {
				temp_pop.append(pop[sub.second[i]]);
			}
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				objs.emplace_back(&temp_pop[i].objective());
			}
			std::vector<int> rank;
			ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				temp_pop[i].setFitness(rank[i]);
			}
			Population<Solution<>> front_pop;
			std::vector<size_t> add_front_inx;
			std::vector<size_t> add_behind_inx;
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				if (temp_pop[i].fitness() == 0) {
					front_pop.append(temp_pop[i]);
					add_front_inx.push_back(i);
				}
				else {
					add_behind_inx.push_back(i);
				}
			}
			//子空间前沿及索引更新，子空间代表解更新
			if (m_subspace_info[idx]->m_subspace_front_sol.empty()) {
				for (size_t i = 0; i < front_pop.size(); ++i) {
					//m_subspace_info[idx]->m_subspace_front_sol.emplace_back(new Solution<>(front_pop[i]));
					m_subspace_info[idx]->m_subspace_front_inx.emplace_back(m_subspace_info[idx]->m_history_inds.size());
					m_subspace_info[idx]->m_history_inds.emplace_back(new Solution<>(front_pop[i]));
					m_subspace_info[idx]->m_sub_freq++;
				}
			}
			else {
				auto temp_front_sols = m_subspace_info[idx]->m_subspace_front_sol;
				auto temp_front_inx = m_subspace_info[idx]->m_subspace_front_inx;
				m_subspace_info[idx]->m_subspace_front_sol.clear();
				m_subspace_info[idx]->m_subspace_front_inx.clear();
				std::vector<size_t> final_add_front_inx;
				std::vector<size_t> temp_add_inx;
				for (size_t i = 0; i < temp_front_sols.size(); ++i) {
					temp_add_inx.push_back(i);
				}
				std::vector<size_t> final_add;
				for (size_t i = 0; i < add_front_inx.size(); ++i) {
					bool flag = false;
					auto& obj2 = temp_pop[add_front_inx[i]].objective();
					for (size_t j = 0; j < temp_front_sols.size(); ++j) {
						auto& obj1 = temp_front_sols[j]->objective();
						auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
						if (ship==Dominance::kDominant) {
							flag = true;
							break;
						}
					}
					if (!flag) {
						final_add = temp_add_inx;
						temp_add_inx.clear();
						final_add_front_inx.push_back(add_front_inx[i]);
						for (size_t j = 0; j < final_add.size(); ++j) {
							if (!temp_pop[add_front_inx[i]].dominate(*temp_front_sols[final_add[j]], pro)) {
								temp_add_inx.push_back(final_add[j]);
							}
						}
					}
				}
				//更新子空间前沿解索引
				for (size_t i = 0; i < final_add_front_inx.size(); ++i) {
					m_subspace_info[idx]->m_subspace_front_inx.push_back(m_subspace_info[idx]->m_history_inds.size());
					m_subspace_info[idx]->m_history_inds.emplace_back(new Solution<>(temp_pop[final_add_front_inx[i]]));
					m_subspace_info[idx]->m_sub_freq++;
				}
				for (size_t i = 0; i < temp_add_inx.size(); ++i) {
					m_subspace_info[idx]->m_subspace_front_inx.push_back(temp_front_inx[temp_add_inx[i]]);
				}
				
			}
			for (size_t i = 0; i < add_behind_inx.size(); ++i) {
				m_subspace_info[idx]->m_history_inds.emplace_back(new Solution<>(temp_pop[add_behind_inx[i]]));
				m_subspace_info[idx]->m_sub_freq++;
			}
			//更新子空间前沿解
			for (size_t i = 0; i < m_subspace_info[idx]->m_subspace_front_inx.size(); ++i) {
				m_subspace_info[idx]->m_subspace_front_sol.emplace_back(m_subspace_info[idx]->m_history_inds[m_subspace_info[idx]->m_subspace_front_inx[i]]);
				m_subspace_info[idx]->m_history_inds[m_subspace_info[idx]->m_subspace_front_inx[i]]->setType(0);
			}

			//更新子空间代表解，如何选代表解，代表解的数量与子空间大小相关
			size_t num_var = CAST_CONOP(pro)->numberVariables();
			//select_num = (num_var + 1) * std::floor(std::pow(subspaceTree().getBoxVolume(idx) / min_volume, 1. / 3));
			size_t total_num = m_subspace_info[idx]->m_subspace_front_sol.size();
			size_t select_num = total_num;
			select_num = (size_t)std::ceil(m_subspace_info[idx]->m_represent_num * std::pow(num_var, 2. / 4));
			m_subspace_info[idx]->m_represent_sol.clear();
			if ( total_num <= select_num) {
				m_subspace_info[idx]->m_subspace_represent_inx = m_subspace_info[idx]->m_subspace_front_inx;
			}
			else {
				bool random_select = true;
				if (random_select) {
					for (size_t j = 0; j < select_num; ++j) {
						size_t inx = (size_t)std::floor(m_subspace_info[idx]->m_subspace_front_sol.size() * rnd->uniform.next());
						m_subspace_info[idx]->m_subspace_represent_inx.push_back(m_subspace_info[idx]->m_subspace_front_inx[inx]);
						//m_subspace_info[idx]->m_represent_sol.emplace_back(m_subspace_info[idx]->m_subspace_front_sol[inx]);
					}
				}
				else {
					std::vector<size_t> selected_inx;
					//先加入某一目标下的最值
					Real min_v = INT16_MAX;
					size_t add_inx=0;
					for (size_t i = 0; i<total_num; ++i) {
						if (m_subspace_info[idx]->m_subspace_front_sol[i]->objective()[0] < min_v) {
							min_v = m_subspace_info[idx]->m_subspace_front_sol[i]->objective()[0];
							add_inx = i;
						}
					}
					selected_inx.push_back(add_inx);
					//先得到未选点索引
					std::vector<size_t> no_select_inx;//未选点索引
					for (size_t i = 0; i < total_num; ++i) {
						if (std::find(selected_inx.begin(), selected_inx.end(), i) == selected_inx.end()) {
							no_select_inx.push_back(i);
						}
					}
					while (selected_inx.size() < select_num) {//选出与当前已选点最远的点
						//计算未选点至已选点的最小距离的最大值
						std::vector<Real> min_dist;
						for (size_t i = 0; i < no_select_inx.size(); ++i) {
							std::vector<Real> temp_dist;
							for (size_t j = 0; j < selected_inx.size(); ++j) {
								auto p1 = m_subspace_info[idx]->m_subspace_front_sol[no_select_inx[i]]->objective();
								auto p2= m_subspace_info[idx]->m_subspace_front_sol[selected_inx[j]]->objective();
								temp_dist.push_back(euclideanDistance(p1.begin(),p1.end(),p2.begin()));
							}
							Real min_d = *std::min_element(temp_dist.begin(),temp_dist.end());
							min_dist.push_back(min_d);
						}
						auto s_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(),min_dist.end()));
						selected_inx.push_back(no_select_inx[s_inx]);
						no_select_inx.erase(no_select_inx.begin()+s_inx);
					}
					for (size_t i = 0; i < selected_inx.size(); ++i) {
						m_subspace_info[idx]->m_subspace_represent_inx.push_back(m_subspace_info[idx]->m_subspace_front_inx[selected_inx[i]]);
					}
				}
			}
			//更新子空间代表解
			for (size_t i = 0; i < m_subspace_info[idx]->m_subspace_represent_inx.size(); ++i) {
				m_subspace_info[idx]->m_represent_sol.emplace_back(m_subspace_info[idx]->m_history_inds[m_subspace_info[idx]->m_subspace_represent_inx[i]]);
			}
			m_subspace_info[idx]->add_flag = "yes";
		}
	}

	void MO_HLC::updateSubSpaceInfo(std::vector<std::shared_ptr<Solution<>>>& pop, Problem *pro, Random *rnd) {
		for (size_t i = 0; i < pop.size(); ++i) {
			auto var = pop[i]->variable().vect();
			auto idx = subspaceTree().getRegionIdx(var);
			//子空间历史解更新
			m_subspace_info[idx]->m_history_inds.emplace_back(pop[i]);
			m_subspace_info[idx]->m_sub_freq++;
			//子空间前沿更新，子空间代表解更新
			if (m_subspace_info[idx]->m_subspace_front_sol.empty()) {
				m_subspace_info[idx]->m_subspace_front_sol.emplace_back(pop[i]);
				m_subspace_info[idx]->m_represent_sol.emplace_back(pop[i]);
			}
			else {
				bool flag = false;
				for (size_t j = 0; j < m_subspace_info[idx]->m_subspace_front_sol.size(); ++j) {
					if (m_subspace_info[idx]->m_subspace_front_sol[j]->dominate(*pop[i], pro->optimizeMode())) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					//更新效率太低，需要改进
					auto temp_front_sols = m_subspace_info[idx]->m_subspace_front_sol;
					m_subspace_info[idx]->m_subspace_front_sol.clear();
					m_subspace_info[idx]->m_subspace_front_sol.emplace_back(pop[i]);
					for (size_t j = 0; j < temp_front_sols.size(); ++j) {
						if (!pop[i]->dominate(*temp_front_sols[j], pro->optimizeMode())) {
							m_subspace_info[idx]->m_subspace_front_sol.emplace_back(temp_front_sols[j]);
						}
					}
					//update represent ind
					if (m_subspace_info[idx]->m_subspace_front_sol.size() <= m_subspace_info[idx]->m_represent_num) {
						m_subspace_info[idx]->m_represent_sol.clear();
						for (size_t j = 0; j < m_subspace_info[idx]->m_subspace_front_sol.size(); ++j) {
							m_subspace_info[idx]->m_represent_sol.emplace_back(m_subspace_info[idx]->m_subspace_front_sol[j]);
						}
					}
					else {
						m_subspace_info[idx]->m_represent_sol.clear();
						for (size_t j = 0; j < m_subspace_info[idx]->m_represent_num; ++j) {
							size_t inx = (size_t)std::floor(m_subspace_info[idx]->m_subspace_front_sol.size() * rnd->uniform.next());
							m_subspace_info[idx]->m_represent_sol.emplace_back(m_subspace_info[idx]->m_subspace_front_sol[inx]);
						}
					}
				}
			}
		}
	}

	void MO_HLC::updateSubspaceInfo(Solution<>& sol) {

	}

	void MO_HLC::updateBasinInfo(const std::vector<std::vector<size_t>>& explore_basin, const std::vector<std::vector<size_t>>& exploit_basin,bool b) {
		m_basin_info.clear();
		m_clusters_MOP.clear();
		////更新cluster
		//for (size_t i = 0; i < explore_basin.size(); ++i) {
		//	m_clusters_MOP.push_back(explore_basin[i]);
		//}
		//for (size_t i = 0; i < exploit_basin.size(); ++i) {
		//	m_clusters_MOP.push_back(exploit_basin[i]);
		//}
		////更新basin
		//for (size_t i = 0; i < m_clusters_MOP.size(); ++i) {
		//	m_basin_info.emplace_back(new BasinInfo_MOP);
		//	//吸引域标签
		//	if (i >= explore_basin.size()) {
		//		m_basin_info.back()->flag = "exploit";
		//	}
		//	//吸引域的探索潜力值
		//	Real explore_potential = 0.;//采用最值或均值
		//	//吸引域的开发潜力值
		//	Real exploit_potential = 0.;//采用最值或均值
		//	if (b) {//取均值
		//		for (size_t j = 0; j < m_clusters_MOP[i].size(); ++j) {
		//			auto tmp2 = getSubspaceInfo(m_clusters_MOP[i][j]).m_explore_potential;
		//			auto tmp3 = getSubspaceInfo(m_clusters_MOP[i][j]).m_exploit_potential;
		//			explore_potential += tmp2;
		//			exploit_potential += tmp3;
		//		}
		//		m_basin_info.back()->m_basin_explore_potential = explore_potential / m_clusters_MOP[i].size();
		//		m_basin_info.back()->m_basin_exploit_potential = exploit_potential / m_clusters_MOP[i].size();
		//		m_basin_info.back()->m_basin_potential = (exploit_potential + explore_potential) / m_clusters_MOP[i].size();
		//	}
		//	else {//取最值
		//		for (size_t j = 0; j < m_clusters_MOP[i].size(); ++j) {
		//			auto tmp2 = getSubspaceInfo(m_clusters_MOP[i][j]).m_explore_potential;
		//			auto tmp3 = getSubspaceInfo(m_clusters_MOP[i][j]).m_exploit_potential;
		//			if (tmp2 >= explore_potential) {
		//				explore_potential = tmp2;
		//			}
		//			if (tmp3 >= exploit_potential) {
		//				exploit_potential = tmp3;
		//			}
		//		}
		//		m_basin_info.back()->m_basin_explore_potential = explore_potential;
		//		m_basin_info.back()->m_basin_exploit_potential = exploit_potential;
		//		m_basin_info.back()->m_basin_potential = exploit_potential + explore_potential;
		//	}
		//	//吸引域内当前最好排序值、当前个体、搜索频率、历史个体
		//	size_t rank = INT16_MAX;
		//	std::vector<std::shared_ptr<Solution<>>> ind;
		//	std::vector<std::shared_ptr<Solution<>>> his_ind;
		//	size_t fre = 0;
		//	for (size_t j = 0; j < m_clusters_MOP[i].size(); ++j) {
		//		auto tmp1 = getSubspaceInfo(m_clusters_MOP[i][j]).m_best_rank;
		//		if (rank > tmp1) {
		//			rank = tmp1;
		//		}
		//		auto sol = getSubspaceInfo(m_clusters_MOP[i][j]).m_curr_sols;
		//		auto his_sol = getSubspaceInfo(m_clusters_MOP[i][j]).m_history_inds;
		//		for (size_t k = 0; k < sol.size(); ++k) {
		//			ind.emplace_back(sol[k]);
		//		}
		//		for (size_t k = 0; k < his_sol.size(); ++k) {
		//			his_ind.emplace_back(his_sol[k]);
		//		}
		//		auto tmp2 = getSubspaceInfo(m_clusters_MOP[i][j]).m_sub_freq;
		//		fre += tmp2;
		//	}
		//	m_basin_info.back()->m_best_rank = rank;
		//	m_basin_info.back()->m_current_indi = ind;
		//	m_basin_info.back()->m_history_inds = his_ind;
		//	m_basin_info.back()->m_basin_freq = fre;
		//	
		//	//吸引域的索引及子空间所属吸引域和类型
		//	m_basin_info.back()->m_cluster_inx = m_basin_info.size() - 1;
		//	for (size_t j = 0; j < m_clusters_MOP[i].size(); ++j) {
		//		m_subspace_info[m_clusters_MOP[i][j]]->idx_cluster.push_back(m_basin_info.size() - 1);
		//		m_subspace_info[m_clusters_MOP[i][j]]->search_flag = m_basin_info.back()->flag;
		//	}
		//	//吸引域包含的子空间
		//	m_basin_info.back()->m_subspace_set = m_clusters_MOP[i];
		//	//吸引域的采样系数
		//	m_basin_info.back()->m_coff_sample = 1;
		//	//吸引域的反馈系数
		//	m_basin_info.back()->m_coff_feedback = 1;
		//}
	}

	//Real MO_HLC::calExploitSpacePotential(size_t idx, std::vector<size_t>& spaces) {
	//	//根据子空间的采样反馈系数、最佳排序值、子目标的最好排序值,以及目标空间的特征影响的潜力值
	//	Real potential;
	//	auto spaceinfo = getSubspaceInfo(idx);
	//	auto sample_feedback = spaceinfo.m_coff_feedback;
	//	//size_t fre = spaceinfo.m_sub_freq;
	//	size_t nd_rank = spaceinfo.m_best_rank;
	//	auto subobjrank = spaceinfo.m_subObj_optima_rank;//目标值越好，rank=1
	//	size_t best_rank = INT16_MAX;
	//	for (size_t i = 0; i < subobjrank.size(); ++i) {
	//		if (best_rank > subobjrank[i]) {
	//			best_rank = subobjrank[i];
	//		}
	//	}
	//	//考察排序值和最好子目标排序
	//	//potential = (1. / (1. + nd_rank)) + (1. / best_rank);
	//	//考察采样反馈、排序值和最好子目标排序
	//	potential = sample_feedback * (1. / (1. + nd_rank)) + (1. / best_rank);
	//	if (std::find(spaces.begin(), spaces.end(), idx) != spaces.end()) {
	//		potential += 0.5;
	//	}

	//	return potential;
	//}

	//Real MO_HLC::calExploreSpacePotential(size_t idx) {
	//	//根据子空间的采样频率，采样反馈系数，子目标的排序值的平均等来确定探索潜力值
	//	//若子空间没有个体，设其子目标值为当前最差值
	//	Real potential;
	//	auto spaceinfo = getSubspaceInfo(idx);
	//	size_t fre = spaceinfo.m_sub_freq;
	//	auto subobjrank = spaceinfo.m_subObj_optima_rank;//目标值越好，rank=1
	//	size_t best_rank = INT16_MAX;
	//	for (size_t i = 0; i < subobjrank.size(); ++i) {
	//		if (best_rank > subobjrank[i]) {
	//			best_rank = subobjrank[i];
	//		}
	//	}
	//	auto sample_feedback = spaceinfo.m_coff_feedback;

	//	//考察搜索频率和最好子目标排序
	//	potential = (1. / (1. + fre)) + (1./best_rank);

	//	//考察采样反馈、搜索频率和最好子目标排序
	//	potential = sample_feedback * (1. / (1. + fre)) + (1. / best_rank);

	//	/*if (std::find(spaces.begin(), spaces.end(), idx) != spaces.end()) {
	//		potential += 0.5;
	//	}*/

	//	return potential;
	//}

	void MO_HLC::spaceDivide(int idx,size_t sub_num) {//to do
		auto bound=subspaceTree().getBox(idx);
		nanoflann::KDTreeSpace<Real> subtree;
		subtree.setInitBox(bound);
		subtree.inputRatioData(std::vector<Real>(sub_num, 1.0 / sub_num));
		subtree.buildIndex();
		std::map<size_t, size_t> add_id;
		subspaceTree().addSubtree(idx, subtree, add_id);
		size_t pre_space_num = numSubspace();
		////先得到划分子空间的历史解和当前解
		//SubspaceInfo_MOP temp_subspace_info = getSubspaceInfo(idx);
		//再增加子空间结构
		for (size_t i = 0; i < sub_num-1; ++i) {
			m_subspace_info.emplace_back(new SubspaceInfo_MOP);
		}
		for (size_t i = pre_space_num; i < m_subspace_info.size(); ++i) {
			m_subspace_info[i]->subspace_inx = i;
			m_subspace_info[i]->idx_cluster = m_subspace_info[idx]->idx_cluster;
			m_subspace_info[i]->m_subspace_granularity = m_subspace_info[idx]->m_subspace_granularity;
			m_subspace_info[pre_space_num]->m_best_rank = m_subspace_info[idx]->m_best_rank;
			//subspaceTree().findNeighbor(i, m_subspace_info[i]->m_sub_neighbors);
		}
		//清空第idx的子空间信息
		m_subspace_info[idx]->m_best_rank = m_subspace_info[idx]->m_best_rank;
		m_subspace_info[idx]->m_linear_neigh_space_inx.clear();
		m_subspace_info[idx]->m_history_inds.clear();
		m_subspace_info[idx]->m_curr_sols.clear();
		m_subspace_info[idx]->m_subspace_front_sol.clear();
		m_subspace_info[idx]->m_subspace_front_inx.clear();
		m_subspace_info[idx]->m_subspace_represent_inx.clear();
		m_subspace_info[idx]->m_front_sol_in_subspace.clear();
		m_subspace_info[idx]->m_history_sub_obj_sol.clear();
		m_subspace_info[idx]->m_gen_front_sols.clear();
		m_subspace_info[idx]->m_sub_neighbors.clear();
		m_subspace_info[idx]->m_represent_sol.clear();
		m_subspace_info[idx]->m_sub_freq = 0;
		//subspaceTree().findNeighbor(idx, m_subspace_info[idx]->m_sub_neighbors);
	    //在算法中使用历史解更新细分子空间信息
	}

	void MO_HLC::spaceSplit(int idx, int dim, Real pos){//to do
		//auto bound = subspaceTree().getBox(idx);
		subspaceTree().splitRegion(idx,&dim,&pos);
		size_t pre_space_num = numSubspace();
		//再增加子空间结构
		m_subspace_info.emplace_back(new SubspaceInfo_MOP);
		m_subspace_info[pre_space_num]->subspace_inx = pre_space_num;
		m_subspace_info[pre_space_num]->idx_cluster = m_subspace_info[idx]->idx_cluster;
		m_subspace_info[pre_space_num]->m_subspace_granularity = m_subspace_info[idx]->m_subspace_granularity;
		m_subspace_info[pre_space_num]->m_best_rank = m_subspace_info[idx]->m_best_rank;
		//清空第idx的子空间信息
		m_subspace_info[idx]->m_best_rank= m_subspace_info[idx]->m_best_rank;
		m_subspace_info[idx]->m_linear_neigh_space_inx.clear();
		m_subspace_info[idx]->m_history_inds.clear();
		m_subspace_info[idx]->m_curr_sols.clear();
		m_subspace_info[idx]->m_subspace_front_sol.clear();
		m_subspace_info[idx]->m_subspace_front_inx.clear();
		m_subspace_info[idx]->m_subspace_represent_inx.clear();
		m_subspace_info[idx]->m_front_sol_in_subspace.clear();
		m_subspace_info[idx]->m_history_sub_obj_sol.clear();
		m_subspace_info[idx]->m_gen_front_sols.clear();
		m_subspace_info[idx]->m_sub_neighbors.clear();
		m_subspace_info[idx]->m_represent_sol.clear();
		m_subspace_info[idx]->m_sub_freq = 0;
		//subspaceTree().findNeighbor(idx, m_subspace_info[idx]->m_sub_neighbors);
		//在算法中使用历史解更新细分子空间信息
	}
	
	void MO_HLC::rankSubspaceByPotential(const std::vector<size_t>& subspaces) {
		

	}

	void MO_HLC::converged(const Solution<>& center) {

	}
}
