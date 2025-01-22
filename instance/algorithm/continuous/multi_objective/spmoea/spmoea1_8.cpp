#include "spmoea1_8.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_8::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_8::run_() {
		initPop(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		while (!terminating()) {
			evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void SPMOEA1_8::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_8::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(getPop().size() + 1);//最后一个为历史解
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j]);
				}
				/*for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j){
					  m_solution[1].push_back(&getPop()[i].getOffspring()[j]);
				}*/
			}
			auto& his_sols = getHisSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(his_sols[i].get());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA1_8::initiVarSpace(Problem *pro) {
		/*************************************/
		/*      MO_HLC搜索空间划分初始化     */
		/*************************************/
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			m_var_boundary.emplace_back(CAST_CONOP(pro)->range(i));
		}
		getMO_HLC().initialVarSpace(m_var_boundary, numVarRegion());
		m_pop_var_range = m_var_boundary;
		m_num_e_e.resize(2);

	}

	void SPMOEA1_8::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
		size_t pop_num = getPopsize();
		Population<Solution<>> new_pop;
		//初始子种群
		SPMOEA_pop temp_pop(pop_num, pro);
		temp_pop.initialize(pro, rnd);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t j = 0; j < pop_num; ++j) {
			temp_pop[j].evaluate(pro, alg);
		}
		SPMOEA::NDSort(temp_pop);
		//初始化子代
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			new_pop.append(temp_pop[j]);
		}
		temp_pop.setPopState("exploit");
		//初始化子种群的历史和每代前沿解
		Population<Solution<>> front_pop;
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			if (temp_pop[j].fitness() == 0) {
				front_pop.append(temp_pop[j]);
				temp_pop.getPopHisFrontSols().emplace_back(std::make_shared<Solution<>>(temp_pop[j]));
			}
		}
		temp_pop.getPopGenFrontSols().emplace_back(std::make_shared<Population<Solution<>>>(front_pop));
		getPop().append(temp_pop);
		getPop().back().setSearchBox(bound);

		//更新历史信息
		SPMOEA::NDSort(new_pop);
		updateHistoryInfo(new_pop, pro);
		SPMOEA::updateObjRange(new_pop, pro);
		SPMOEA::updateNewPop(new_pop);
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		SPMOEA::updateSubspaceFrontSol(new_pop, pro, rnd);
		updateFrontSpace();
		setFrontLastGens(1);
		updateEE(pop_num, 0);
		//根据子空间划分聚类
		clusterSubspace();
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
	}

	int SPMOEA1_8::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		/********************************************************************************
							            子种群资源分配
		********************************************************************************/
		//根据子种群占有的前沿子空间数、种群个体数决定，最后一个子种群采用初始种群数
		std::vector<size_t> assign_pop_resource;
		PopResourceAssign(assign_pop_resource, pro);
		/********************************************************************************
											 子种群演化
		********************************************************************************/
		//1.子代中探索与开发的比例；2.子代中探索与开发的方式
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			bool m_assign_resource = true;//子种群是否进行资源分配
			if (!m_assign_resource) {
				for (auto& i : assign_pop_resource) {
					i = getPopsize();
				}
			}
			std::vector<int> interactive_type;
			auto front_link_spaces = getFrontRegionLinkSpace();
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (i < front_link_spaces.size()) {
					interactive_type.push_back(1);
				}
				else {
					interactive_type.push_back(1);
				}
			}
			//根据联通情况采样
			generateOffspring(pro, rnd, assign_pop_resource, interactive_type);
		}
		else {
			//generating solutions by learning methods
			//multi-popution and prediction
			//代内排序预测，代际预测更新
			//子种群演化历史轨迹
		}
		Population<Solution<>> offspring_pop;
		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); j++) {
				tag = getPop()[i].getOffspring()[j].evaluate(pro, alg);
				if (tag != EvaluationTag::kNormalEval)
					break;
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}
		
		SPMOEA::NDSort(offspring_pop);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		updateSubspaceFrontSol(offspring_pop, pro, rnd);//使用子代更新子空间前沿解
		
		auto pre_front_space = getFrontSpace();
		updateFrontSpace();//使用子种群子代更新子空间信息
		//auto middle_front_space= getFrontSpace();
		////看最后一个子种群是否找到新的前沿子空间
		//SPMOEA::NDSort(last_offspring_pop);
		//SPMOEA::updateHistoryInfo(last_offspring_pop, pro);
		//updateSubspaceFrontSol(last_offspring_pop, pro, rnd);//使用子代更新子空间前沿解
		//updateFrontSpace();
		//auto final_front_space = getFrontSpace();
		/*std::vector<size_t> explore_new_front_space;
		for (size_t i = 0; i < final_front_space.size(); ++i) {
			if (std::find(middle_front_space.begin(), middle_front_space.end(), final_front_space[i]) == middle_front_space.end()) {
				explore_new_front_space.push_back(final_front_space[i]);
			}
		}*/
		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		//检测前沿子空间是否持续一定代数未变
		auto flag = ifFrontChanged(pre_front_space);
		/*if (!flag && getFrontLastGens() >= 3) {
			spaceSubdivision(pro, rnd);
			updateFrontSpace();
			setFrontLastGens(1);
			clusterSubspace();
		}*/
		if (m_divide_iteration % 2 == 0) {
			spaceSubdivision(pro, rnd);
			updateFrontSpace();
			setFrontLastGens(1);
			clusterSubspace();
		}

		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		/*Population<Solution<>> offspring_pop;
		for (size_t i = 0; i < pre_offspring_pop.size(); ++i) {
			offspring_pop.append(pre_offspring_pop[i]);
		}
		for (size_t i = 0; i < last_offspring_pop.size(); ++i) {
			offspring_pop.append(last_offspring_pop[i]);
		}*/
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				auto v = getMO_HLC().subspaceTree().getBoxVolume(i);
				front_ratio += (v / total_volume);
				front_count++;
			}
		}
		setFrontSpaceRatio(front_ratio);
		/********************************************************************************************
									子种群淘汰选择、合并、生成
		*********************************************************************************************/
		//让子种群的演化具有连续性
		//updatePop(pre_offspring_pop,last_offspring_pop,pro,alg,rnd);
		assignPops(offspring_pop, pro,alg, rnd);
		/*for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			getPop()[i].envirSelection(pro, rnd, 0);
			getPop()[i].iteration()++;
		}*/
		//更新最后一个探索子种群：在剩余的子空间中初始化种群
		//当最后一个子种群的前沿在已有前沿子空间时，重新初始化
		//在分离前沿子空间连通集的中间区域初始化新的子种群

		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA1_8::assignPops(Population<Solution<>>& off_pop, Problem *pro,Algorithm *alg,Random *rnd) {
		//先对前沿子空间聚类
		auto front_spaces = getFrontSpace();
		auto front_clusters = clusterFrontSpace(front_spaces);
		getFrontRegionLinkSpace() = front_clusters;
		//找出前沿子空间连通集及其一级邻域以及非邻域
		std::vector<std::vector<size_t>> temp_clusters;
		std::vector<size_t> all_neighs;
		for (size_t i = 0; i < front_clusters.size(); ++i) {
			auto temp = front_clusters[i];
			std::vector<size_t> neighs = temp;
			for (size_t j = 0; j < temp.size(); ++j) {
				auto nei = getMO_HLC().getSubspaceInfo(temp[j]).m_sub_neighbors;
				for (auto jj : nei) {
					if (std::find(neighs.begin(), neighs.end(), jj) == neighs.end()) {
						neighs.push_back(jj);
					}
				}
			}
			/*for (auto jj : neighs) {
				all_neighs.push_back(jj);
			}
			temp_clusters.emplace_back(neighs);*/
			//只考虑前沿子空间
			for (auto jj : temp) {
				all_neighs.push_back(jj);
			}
			temp_clusters.emplace_back(temp);
		}
		std::vector<size_t> other_spaces;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (std::find(all_neighs.begin(), all_neighs.end(), i) == all_neighs.end()) {
				other_spaces.push_back(i);
			}
		}
		//找出父代个体依附于哪一个连通集
		Population<Solution<>> parent_pop;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); j++) {
				parent_pop.append(getPop()[i][j]);
			}
		}
		std::vector<size_t> parent_attach_cluster(parent_pop.size(), temp_clusters.size());
		for (size_t i = 0; i < parent_pop.size(); ++i) {
			auto sol = parent_pop[i].variable().vect();
			auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			for (size_t j = 0; j < temp_clusters.size(); ++j) {
				if (std::find(temp_clusters[j].begin(), temp_clusters[j].end(), inx) != temp_clusters[j].end()) {
					parent_attach_cluster[i] = j;
					break;
				}
			}
		}
		//再找出子代个体依附于哪一个连通集
		std::vector<size_t> off_attach_cluster(off_pop.size(), temp_clusters.size());
		for (size_t i = 0; i < off_pop.size(); ++i) {
			auto sol = off_pop[i].variable().vect();
			auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			for (size_t j = 0; j < temp_clusters.size(); ++j) {
				if (std::find(temp_clusters[j].begin(), temp_clusters[j].end(), inx) != temp_clusters[j].end()) {
					off_attach_cluster[i] = j;
					break;
				}
			}
		}
		//基于每一个前沿连通集进行种群更新
		auto pop_size = getPopsize();
		auto& bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < getPop().size(); ++i) {
			while (getPop()[i].size() > 0) {
				getPop()[i].remove(getPop()[i].end() - 1);
			}
			while (getPop()[i].getOffspring().size() > 0) {
				getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
			}
		}
		while (getPop().size() < temp_clusters.size()) {
			SPMOEA_pop new_pop(0, pro);
			new_pop.setSearchBox(bound);
			getPop().append(new_pop);
		}
		while (getPop().size() > temp_clusters.size()) {
			getPop().remove(getPop().end()-1);
		}
		for (size_t i = 0; i < getPop().size(); ++i) {
			getPop()[i].setPopState("exploit");
		}
		if (other_spaces.size()>0) {
			//增加探索种群
			SPMOEA_pop new_pop(0, pro);
			new_pop.setPopState("explore");
			new_pop.setSearchBox(bound);
			getPop().append(new_pop);
			temp_clusters.emplace_back(other_spaces);
		}
		for (size_t i = 0; i < parent_attach_cluster.size(); ++i) {
			auto cluster_inx = parent_attach_cluster[i];
			getPop()[cluster_inx].append(parent_pop[i]);
		}
		for (size_t i = 0; i < off_attach_cluster.size(); ++i) {
			auto cluster_inx = off_attach_cluster[i];
			getPop()[cluster_inx].append(off_pop[i]);
		}
		//种群更新，环境选择
		size_t M = CAST_CONOP(pro)->numberObjectives();
		size_t con = CAST_CONOP(pro)->numberConstraints();
		size_t var = CAST_CONOP(pro)->numberVariables();
		Population<Solution<>> new_ind;
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState()=="explore") {//探索子种群
				if (front_clusters.size() == 1) {
					//auto manifold_dist = calSpaceManifoldDist(front_clusters[0]);
					////找出相距最远的两个子空间
					//size_t inx1, inx2;
					//std::vector<size_t> temp_max;
					//for (size_t j = 0; j < manifold_dist.size(); ++j) {
					//	temp_max.push_back(*std::max_element(manifold_dist[j].begin(), manifold_dist[j].end()));
					//}
					//auto index1 = std::distance(temp_max.begin(), std::max_element(temp_max.begin(), temp_max.end()));
					//auto index2 = std::distance(manifold_dist[index1].begin(), std::max_element(manifold_dist[index1].begin(), manifold_dist[index1].end()));
					//inx1 = front_clusters[0][index1];
					//inx2 = front_clusters[0][index2];
					////找出位于两个子空间邻域且为非前沿子空间的子空间

					////将位于这些子空间的个体当做探索子种群个体


				}
				else {
					//找出连通子空间最近的子空间
					auto select_pair_spaces = findCloseSpaces(front_clusters, rnd);
					Real min_dist = INT16_MAX;
					size_t inx = 0;
					for (size_t j = 0; j < select_pair_spaces.size(); ++j) {
						if (std::get<2>(select_pair_spaces[j]) < min_dist) {
							min_dist = std::get<2>(select_pair_spaces[j]);
							inx = j;
						}
					}
					//找出位于两个子空间中心之间的个体当做更新的探索种群
					auto box1 = getMO_HLC().subspaceTree().getBox(std::get<0>(select_pair_spaces[inx]));
					auto box2 = getMO_HLC().subspaceTree().getBox(std::get<1>(select_pair_spaces[inx]));
					std::vector<Real> center1,center2,center3;
					//std::vector<std::pair<Real, Real>> initial_bound(box1.size(),std::make_pair<>(0.,0.));
					for (size_t j = 0; j < box1.size(); ++j) {
						center1.push_back((box1[j].first+box1[j].second)/2);
						center2.push_back((box2[j].first + box2[j].second) / 2);
						center3.push_back((center1.back()+center2.back())/2);
						/*initial_bound[j].first = center1.back() < center2.back() ? center1.back() : center2.back();
						initial_bound[j].second = center1.back() > center2.back() ? center1.back() : center2.back();
						if (initial_bound[j].first == initial_bound[j].second) {
							initial_bound[j].first = box1[j].first < box2[j].first ? box1[j].first : box2[j].first;
							initial_bound[j].second = box1[j].second > box2[j].second ? box1[j].second : box2[j].second;
						}*/
						/*initial_bound[j].first = box1[j].first < box2[j].first ? box1[j].first : box2[j].first;
						initial_bound[j].second = box1[j].second > box2[j].second ? box1[j].second : box2[j].second;*/
					}
					size_t space_inx = getMO_HLC().subspaceTree().getRegionIdx(center3);
				    //找出邻域
					std::vector<size_t> neigh_space;
					auto neis = getMO_HLC().getSubspaceInfo(space_inx).m_sub_neighbors;
					neis.push_back(space_inx);
					for (auto jj : neis) {
						if (std::find(all_neighs.begin(), all_neighs.end(), jj) == all_neighs.end()) {
							neigh_space.push_back(jj);
						}
					}
					std::vector<size_t> select_flag(getPop()[i].size(),0);
					for (size_t j = 0; j < getPop()[i].size(); ++j) {
						auto sol = getPop()[i][j].variable().vect();
						auto sp_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
						if (std::find(neigh_space.begin(), neigh_space.end(), sp_inx) != neigh_space.end()) {
							select_flag[j] = 1;
						}
					}
					for (int j = select_flag.size()-1; j >=0; --j) {
						if (select_flag[j] == 0) {
							getPop()[i].remove(getPop()[i].begin()+j);
						}
					}
					while (getPop()[i].size() < pop_size) {
						std::vector<Real> new_sol;
						if (neigh_space.size() > 0) {
							auto index = (size_t)std::floor(neigh_space.size() * rnd->uniform.next());
							auto box = getMO_HLC().subspaceTree().getBox(neigh_space[index]);
							for (size_t k = 0; k < box.size(); ++k) {
								new_sol.push_back(box[k].first + (box[k].second - box[k].first) * rnd->uniform.next());
							}
						}
						else {

						}
						Solution<> ind(M, con, var);
						ind.variable().vect() = new_sol;
						ind.evaluate(pro,alg,true);
						getPop()[i].append(ind);
						new_ind.append(ind);
					}
				}
			}
			if(getPop()[i].size() > pop_size) {
				//从中选pop_size个个体
				SPMOEA::NDSort(getPop()[i]);
				Population<Solution<>> first_pop;
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					if (getPop()[i][j].fitness() == 0) {
						first_pop.append(getPop()[i][j]);
					}
				}
				size_t com_rank = 1;
				while (first_pop.size() < pop_size) {
					for (size_t j = 0; j < getPop()[i].size(); ++j) {
						if (getPop()[i][j].fitness() == com_rank) {
							first_pop.append(getPop()[i][j]);
							if (first_pop.size() >= pop_size) {
								break;
							}
						}
					}
					com_rank++;
				}
				std::vector<size_t> sele_inx;
				if (first_pop.size() > pop_size) {
					sele_inx = selectMaxMinFromFront(first_pop,pop_size);
					std::sort(sele_inx.begin(),sele_inx.end());
					bool repeat_ele = false;
					for (size_t j= 0; j < sele_inx.size()-1; ++j) {
						if (sele_inx[j] == sele_inx[j + 1]) {
							repeat_ele = true;
						}
					}
				}
				else {
					for (size_t j = 0; j < first_pop.size(); ++j) {
						sele_inx.push_back(j);
					}
				}
				while (getPop()[i].size() > pop_size) {
					getPop()[i].remove(getPop()[i].end() - 1);
				}
				for (size_t j = 0; j < sele_inx.size(); ++j) {
					getPop()[i][j] = first_pop[sele_inx[j]];
				}
			}
			else if(getPop()[i].size() < pop_size) {//在子空间生成个体
				Population<Solution<>> add_pop;
				auto spaces = temp_clusters[i];
				for (size_t j = 0; j < pop_size - getPop()[i].size(); ++j) {
					auto index = (size_t)std::floor(spaces.size() * rnd->uniform.next());
					auto box = getMO_HLC().subspaceTree().getBox(spaces[index]);
					std::vector<Real> new_sol;
					for (size_t k = 0; k < var; ++k) {
						new_sol.push_back(box[k].first + (box[k].second - box[k].first) * rnd->uniform.next());
					}
					Solution<> ind(M, con, var);
					ind.variable().vect() = new_sol;
					add_pop.append(ind);
				}
				add_pop.evaluate(pro, alg);
				for (size_t j = 0; j < add_pop.size(); ++j) {
					getPop()[i].append(add_pop[j]);
					new_ind.append(add_pop[j]);
				}
			}
			//更新子代
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring().append(getPop()[i][j]);
				getPop()[i].getOffspring().append(getPop()[i][j]);
			}
		}
		if (new_ind.size()>0) {
			SPMOEA::NDSort(new_ind);
			SPMOEA::updateHistoryInfo(new_ind, pro);
			updateSubspaceFrontSol(new_ind, pro, rnd);//使用新解更新子空间前沿解
		}
	}

	void SPMOEA1_8::updatePop(Population<Solution<>>& off_pop1, Population<Solution<>>& off_pop2, Problem *pro, Algorithm *alg, Random *rnd) {
		//让子种群的演化具有连续性，父代不合并
		//子代个体离哪个子种群近，就归属到哪个子种群
		//子代进入到某一个子种群所在子空间的邻域，则对这个子种群进行淘汰选择，
		//剩余子代与探索种群合并演化，若有在前沿子空间的，需要生成新的探索种群
		//先对前沿子空间聚类
		auto front_clusters = clusterFrontSpace(getFrontSpace());
		//再找出种群依附于哪一个连通集
		std::vector<size_t> attach_cluster;
		std::vector<size_t> lost_clusters;//在分离的新的前沿子空间及其邻域生成多种群
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			//父代个体在每个连通集的个数
			std::vector<size_t> num_in_space(front_clusters.size(),0);
			for (size_t j = 0; j < front_clusters.size(); ++j) {
				for (size_t k = 0; k < getPop()[i].size(); ++k) {
					auto sol = getPop()[i][k].variable().vect();
					auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					if (std::find(front_clusters[j].begin(), front_clusters[j].begin(), inx) != front_clusters[j].end()) {
						num_in_space[j]++;
					}
				}
			}
			attach_cluster.push_back(std::distance(num_in_space.begin(),std::max_element(num_in_space.begin(),num_in_space.end())));
		}
		for (size_t i = 0; i < front_clusters.size(); ++i) {
			if (std::find(attach_cluster.begin(), attach_cluster.end(), i) == attach_cluster.end()) {
				lost_clusters.push_back(i);
			}
		}
		//父代种群位于的子空间
		std::vector<std::vector<size_t>> pop_spaces;
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			std::vector<size_t> temp;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(getPop()[i][j].variable().vect());
				if (std::find(temp.begin(),temp.end(),space_inx)==temp.end()) {
					temp.push_back(space_inx);
				}
			}
			pop_spaces.emplace_back(temp);
		}
		//挑选父代
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			std::vector<size_t> temp;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(getPop()[i][j].variable().vect());
				if (std::find(temp.begin(), temp.end(), space_inx) == temp.end()) {
					temp.push_back(space_inx);
				}
			}
			pop_spaces.emplace_back(temp);
		}
		//挑选子代
		

		//生成探索子种群，如果探索子种群进入到已有子种群范围内，则需要重新生成
		if (lost_clusters.size() > 0) {

		}
		//取两个空间集的交集及其邻域作为子种群的演化
		//子代依附的种群
		std::vector<std::vector<size_t>> off_assign(getPop().size() - 1);
		for (size_t i = 0; i < off_pop1.size(); ++i) {
			//子代分别到不同的前沿子连通集的距离，选最小的依附
			//先看是否在某一个子空间内
			auto& sol = off_pop1[i].variable().vect();
			auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			std::vector<size_t> flag(getPop().size() - 1, 0);
			for (size_t j = 0; j < pop_spaces.size(); ++j) {
				if (std::find(pop_spaces[j].begin(), pop_spaces[j].end(), inx) != pop_spaces[j].end()) {
					flag[j] = 1;
					break;
				}
			}
			if (std::find(flag.begin(), flag.end(), 1) == flag.end()) {
				//再看离亲代距离的大小
				std::vector<Real> min_dist;
				for (size_t j = 0; j < getPop().size() - 1; ++j) {
					std::vector<Real> temp_dist;
					for (size_t k = 0; k < getPop()[j].size(); ++k) {
						auto dist = euclideanDistance(sol.begin(), sol.end(), getPop()[j][k].variable().vect().begin());
						temp_dist.push_back(dist);
					}
					min_dist.push_back(*std::min_element(temp_dist.begin(),temp_dist.end()));
				}
				auto pop_inx= std::distance(min_dist.begin(), std::min_element(min_dist.begin(), min_dist.end()));
				off_assign[pop_inx].push_back(i);
			}
			else {
				auto pop_inx = std::distance(flag.begin(),std::find(flag.begin(),flag.end(),1));
				off_assign[pop_inx].push_back(i);
			}
		}

		//根据分配的子代数更新各个子种群子代
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			while (getPop()[i].getOffspring().size() - getPop()[i].size() > off_assign[i].size()) {
				getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
			}
			while (getPop()[i].getOffspring().size() - getPop()[i].size() < off_assign[i].size()) {
				Solution<> temp_ind(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
				getPop()[i].getOffspring().append(temp_ind);
			}
		}
		
		//为子代赋值
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			for (size_t j = 0; j < off_assign[i].size(); ++j) {
				getPop()[i].getOffspring()[j].variable() = off_pop1[off_assign[i][j]].variable();
			}
		}


		auto cur_clusters = getMO_HLC().getClusters();
		//提取原始种群个体
		std::map<size_t, std::vector<std::shared_ptr<Solution<>>>> space_inds;//每个子空间含有的当前个体
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(getPop()[i][j].variable().vect());
				if (space_inds[space_inx].empty()) {
					std::vector<std::shared_ptr<Solution<>>> temp_inds;
					space_inds.insert(std::make_pair(space_inx, temp_inds));
				}
				space_inds[space_inx].emplace_back(std::make_shared<Solution<>>(getPop()[i][j]));
			}
		}
		//先删除之前的种群
		while (getPop().size() > 0) {
			getPop().remove(getPop().end() - 1);
		}
		//再根据聚类的结果，将原有个体分配到新种群
		for (size_t i = 0; i < cur_clusters.size(); ++i) {
			size_t count = 0;//新的类中含有的个体数
			size_t pop_size = getPopsize();
			for (size_t j = 0; j < cur_clusters[i].size(); ++j) {
				for (size_t k = 0; k < cur_clusters[i][j].size(); ++k) {
					count += space_inds[cur_clusters[i][j][k]].size();
				}
			}
			SPMOEA_pop temp_pop(pop_size, pro);
			if (count > pop_size) {//聚类子空间个体充足，从子空间所有个体中选择精英个体
				SPMOEA_pop temp_pop0(count, pro);
				size_t count0 = 0;
				for (size_t j = 0; j < cur_clusters[i].size(); ++j) {
					for (size_t k = 0; k < cur_clusters[i][j].size(); ++k) {
						for (size_t p = 0; p < space_inds[cur_clusters[i][j][k]].size(); ++p) {
							temp_pop0[count0] = *space_inds[cur_clusters[i][j][k]][p];
							count0++;
						}
					}
				}
				SPMOEA::NDSort(temp_pop0);
				//分层加入个体
				auto select_inx = selectIndiFromSpace(temp_pop0, pop_size, pro, rnd);
				for (size_t j = 0; j < select_inx.size(); ++j) {
					temp_pop[j] = temp_pop0[select_inx[j]];
				}
			}
			else {//聚类子空间当前个体数不足时
				size_t cc = 0;//子空间现有当前个体数
				std::vector<size_t> cluster_space_inx;//子空间聚类的所有子空间
				for (size_t j = 0; j < cur_clusters[i].size(); ++j) {
					for (size_t k = 0; k < cur_clusters[i][j].size(); ++k) {
						cluster_space_inx.push_back(cur_clusters[i][j][k]);
						for (size_t p = 0; p < space_inds[cur_clusters[i][j][k]].size(); ++p) {
							temp_pop[cc] = *space_inds[cur_clusters[i][j][k]][p];
							cc++;
						}
					}
				}
				if (count < pop_size) {
					//在子空间内额外采样: 随机 or 精英采样
					SPMOEA_pop temp_pop1(pop_size - count, pro);
					for (size_t j = 0; j < pop_size - count; ++j) {
						size_t inx = std::floor(cluster_space_inx.size() * rnd->uniform.next());
						auto bound = getMO_HLC().subspaceTree().getBox(cluster_space_inx[inx]);
						for (size_t k = 0; k < bound.size(); ++k) {
							temp_pop1[j].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
						}
					}
					temp_pop1.evaluate(pro, alg);
					SPMOEA::NDSort(temp_pop1);
					//更新子空间信息
					updateVarSpaceInfo(temp_pop1, pro, rnd);
					//更新历史前沿
					updateHistoryInfo(temp_pop1, pro);
					for (size_t j = 0; j < temp_pop1.size(); ++j) {
						temp_pop[cc + j] = temp_pop1[j];
					}
				}
			}
			//初始化子代
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring()[j] = temp_pop[j];
				temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			}
			getPop().append(temp_pop);
		}
		//更新最后一个探索子种群：在剩余的子空间中初始化种群
		//当最后一个子种群的前沿在已有前沿子空间时，重新初始化
		//在分离前沿子空间连通集的中间区域初始化新的子种群
		//如果子种群找到新的前沿子空间，则在特定区域生成新的探索子种群，若没有找到，则去除父代及子代在已有前沿的个体，并在特定区域补齐
	}

	void SPMOEA1_8::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			auto pop_state = getPop()[i].getPopState();
			auto search_box = getPop()[i].getPopSearchBox();
			//种群历史前沿解在搜索域的个数
			size_t in_count = 0;
			auto& pop_his_front_sols = getPop()[i].getPopHisFrontSols();
			for (size_t j = 0; j < pop_his_front_sols.size(); ++j) {
				if (indInRange(search_box, pop_his_front_sols[j]->variable().vect())) {
					in_count++;
				}
			}
			Real in_ratio = (Real)in_count / pop_his_front_sols.size();
			front_bound_ratio.push_back(in_ratio);
			if (in_ratio < 0.1 && pop_state == "active") {//参数
				getPop()[i].setPopState("sleep");
			}
		}
	}

	bool SPMOEA1_8::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_8::updateSubPopSpace(size_t pop_inx, Problem *pro, Random *rnd) {
		//根据新种群的信息，更新子种群所在搜索域的子空间排序值
		//更新子空间在子种群中吸引域内的最好rank值
		if (getPop()[pop_inx].getPopState() == "active") {
			Population<Solution<>> temp_off_pop;
			for (size_t j = 0; j < getPop()[pop_inx].getOffspring().size() - getPop()[pop_inx].size(); ++j) {
				temp_off_pop.append(getPop()[pop_inx].getOffspring()[j]);
			}
			getMO_HLC().updateSubspaceInfo(temp_off_pop, pro, rnd);

			std::vector<std::shared_ptr<Solution<>>> temp_pop1;
			for (size_t j = 0; j < getMO_HLC().getCluster(pop_inx)[0].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_represent_sol.size(); ++k) {
					temp_pop1.emplace_back(getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_represent_sol[k]);
				}
			}
			std::vector<std::vector<Real>*> objs;
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				objs.emplace_back(&temp_pop1[j]->objective());
			}
			std::vector<int> rank;
			int num_layer = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				temp_pop1[j]->setFitness(rank[j]);
			}
			//更新种群搜索域子空间最好排序值
			for (size_t j = 0; j < getMO_HLC().getCluster(pop_inx)[0].size(); ++j) {
				getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = INT16_MAX;
			}
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				auto temp_sol = temp_pop1[j]->variable().vect();
				auto inx = getMO_HLC().subspaceTree().getRegionIdx(temp_sol);
				if (temp_pop1[j]->fitness() < getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank) {
					getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank = temp_pop1[j]->fitness();
				}
			}
			//对于没有个体的子空间，找邻域最好的排序值
			std::vector<size_t> cluster_spaces = getMO_HLC().getCluster(pop_inx)[0];
			for (size_t j = 0; j < cluster_spaces.size(); ++j) {
				if (getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_best_rank == INT16_MAX) {
					auto neigh_space = getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_neighbors;
					std::vector<int> neigh_ranks;
					for (auto iter = neigh_space.begin(); iter != neigh_space.end(); ++iter) {
						size_t inx = *iter;
						if (std::find(cluster_spaces.begin(), cluster_spaces.end(), inx) != cluster_spaces.end()) {
							if (getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank < INT16_MAX) {
								neigh_ranks.push_back(getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank);
							}
						}
					}
					if (neigh_ranks.empty()) {
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = num_layer;
					}
					else if (neigh_ranks.size() == 1) {
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
					}
					else {
						int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
						int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = std::ceil((min_rank + max_rank) / 2) + 1;
					}
				}
			}
		}
	}

	void SPMOEA1_8::PopResourceAssign(std::vector<size_t>& assign_pop_resource, Problem *pro) {
		size_t max_pop_size = getPopsize();
		//updateFrontRegionLinkSpace();
		//auto front_link_spaces = getFrontRegionLinkSpace();
		//size_t M = CAST_CONOP(pro)->numberObjectives();
		//std::vector<size_t> front_space_num;
		//for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
		//	auto spaces = getMO_HLC().getCluster(i)[0];
		//	size_t count = 0;
		//	for (size_t j = 0; j < spaces.size(); ++j) {
		//		if (getMO_HLC().getSubspaceInfo(spaces[j]).m_best_rank == 0) {
		//			count++;
		//		}
		//	}
		//	front_space_num.push_back(count);
		//}
		//for (size_t i = 0; i < front_space_num.size(); ++i) {
		//	assign_pop_resource.push_back(front_space_num[i]+5);
		//	//assign_pop_resource.push_back(max_pop_size);
		//}
		//assign_pop_resource.push_back(max_pop_size);
		for (size_t i = 0; i < getPop().size(); ++i) {
			assign_pop_resource.push_back(max_pop_size);
		}
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_8::generateOffspring(Problem *pro, Random *rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
		auto search_bound = CAST_CONOP(pro)->boundary();
		//前面的种群在自己的连通集之间扩展
		size_t total_num = 0;
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t exploit_num = pop_resource[i];
			total_num += exploit_num;
			//子代个体数量处理
			while (getPop()[i].getOffspring().size() - getPop()[i].size() > pop_resource[i]) {
				getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
			}
			while (getPop()[i].getOffspring().size() - getPop()[i].size() < pop_resource[i]) {
				Solution<> temp_ind(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
				getPop()[i].getOffspring().append(temp_ind);
			}
			//在子区域前排个体间开发,开发方式选择与邻域个体之间交互
			if (type[i] == 1) {
				//auto all_off = sampleByGA(getPop()[i], exploit_num, pro, rnd);
				size_t kk = 2;
				auto bound = getPop()[i].getPopSearchBox();
				auto all_off = sampleByDE(getPop()[i], bound, exploit_num, kk, pro, rnd);
				//为子代赋值
				for (size_t k = 0; k < all_off.size(); ++k) {
					getPop()[i].getOffspring()[k].variable().vect() = all_off[k];
					getPop()[i].getOffspring()[k].setCounter(0);
				}
			}
			else if (type[i] == 2) {//前沿子空间连通集交互，连通集中采样
				auto front_link_spaces = getFrontRegionLinkSpace();
				auto space_count = getMO_HLC().numSubspace();
				auto bound = CAST_CONOP(pro)->boundary();
				//邻域前沿子空间交互
				std::vector<std::vector<Real>> all_off;
				auto off = sampleInLinkSpaces(front_link_spaces[i], bound, exploit_num, pro, rnd);
				for (auto& jj : off) {
					all_off.emplace_back(jj);
				}
				for (size_t j = 0; j < all_off.size(); ++j) {
					getPop()[i].getOffspring()[j].variable().vect() = all_off[j];
					getPop()[i].getOffspring()[j].setCounter(0);
				}
				//if (front_link_spaces.size() == 1) {//子空间连通
				//	//邻域前沿子空间交互
				//	std::vector<std::vector<Real>> all_off;
				//	auto off = sampleInLinkSpaces(front_link_spaces[0], bound, getPop()[i].size(), pro, rnd);
				//	for (auto& jj : off) {
				//		all_off.emplace_back(jj);
				//	}
				//	for (size_t j = 0; j < all_off.size(); ++j) {
				//		getPop().back().getOffspring()[j].variable().vect() = all_off[j];
				//		getPop().back().getOffspring()[j].setCounter(0);
				//	}
				//}
				//else {
				//	//前沿子连通集内部交互
				//	std::vector<std::vector<Real>> all_off;
				//	size_t off_count = 0;
				//	for (size_t j = 0; j < front_link_spaces.size(); ++j) {
				//		//子连通集内部交互
				//		auto off2 = sampleInLinkSpaces(front_link_spaces[j], bound, front_link_spaces[j].size(), pro, rnd);
				//		for (size_t k = 0; k < off2.size(); ++k) {
				//			getPop().back().getOffspring()[off_count].variable().vect() = off2[k];
				//			getPop().back().getOffspring()[off_count].setCounter(0);
				//			off_count++;
				//		}
				//	}

				//	//子连通集之间的交互（填补空缺）
				//	auto select_pair_spaces = findCloseSpaces(front_link_spaces, rnd);
				//	//分别在每一个子连通集上进行最小距离子空间的交互
				//	for (size_t j = 0; j < select_pair_spaces.size(); ++j) {
				//		size_t inx1 = std::get<0>(select_pair_spaces[j]);
				//		size_t inx2 = std::get<1>(select_pair_spaces[j]);
				//		int method = 1;
				//		if (method == 1) {//在这两个子空间采用向量方式采样
				//			//子空间前沿个体交互
				//			auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_sub;
				//			auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_sub;
				//			size_t s1 = std::floor(front_sols1.size() * rnd->uniform.next());
				//			size_t s2 = std::floor(front_sols2.size() * rnd->uniform.next());
				//			auto p1 = front_sols1[s1]->variable().vect();
				//			auto p2 = front_sols2[s2]->variable().vect();

				//			for (size_t k = 0; k < 5; ++k) {
				//				auto off = vectorOperator(p1, p2, search_bound, 0, rnd);
				//				getPop().back().getOffspring()[off_count].variable().vect() = off;
				//				getPop().back().getOffspring()[off_count].setCounter(0);
				//				off_count++;
				//			}
				//		}
				//		else if (method == 2) {//在两个子空间算子方式采样
				//			auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_sub;
				//			auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_sub;
				//			//将这些个体组成一个种群，产生新解
				//			SPMOEA_pop temp_pop(0, pro);
				//			for (size_t j = 0; j < front_sols1.size(); ++j) {
				//				temp_pop.append(*front_sols1[j]);
				//			}
				//			for (size_t j = 0; j < front_sols2.size(); ++j) {
				//				temp_pop.append(*front_sols2[j]);
				//			}
				//			for (size_t j = 0; j < temp_pop.size(); ++j) {
				//				temp_pop.getOffspring().append(temp_pop[j]);
				//			}
				//			for (size_t j = 0; j < temp_pop.size(); ++j) {
				//				temp_pop.getOffspring().append(temp_pop[j]);
				//			}
				//			size_t kk = 2;
				//			std::vector<std::vector<Real>> off;
				//			if (temp_pop.size() < 3) {
				//				off = sampleByGA(temp_pop, bound, 5, pro, rnd);
				//			}
				//			else {
				//				off = sampleByDE(temp_pop, bound, 5, 2, pro, rnd);
				//			}
				//			for (size_t k = 0; k < off.size(); ++k) {
				//				getPop().back().getOffspring()[off_count].variable().vect() = off[k];
				//				getPop().back().getOffspring()[off_count].setCounter(0);
				//				off_count++;
				//			}
				//		}
				//	}
				//}
			}
		}
		updateEE(0, total_num);
	}

	void SPMOEA1_8::spaceSubdivision(Problem *pro, Random *rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		auto front_spaces = getFrontSpace();
		//细分子空间
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			size_t min_num = getSplitGranularity();
			// min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
			if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
				//先得到子空间历史解用于更新细分的子空间的信息
				int dim = findSplitDim(sp, pro);
				auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
				Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
				splitSpace(sp, 2, dim, pos, false, pro, rnd);
			}
		}

		//子空间体积占比
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_8::clusterSubspace() {
		getMO_HLC().getClusters().clear();
		auto front_space = getFrontSpace();
		auto front_cluster = clusterFrontSpace(front_space);
		//找出前沿子空间的一级邻域
		std::vector<std::vector<size_t>> temp_cluster;
		for (size_t i = 0; i < front_cluster.size(); ++i) {
			auto temp = front_cluster[i];
			temp_cluster.emplace_back(temp);
			std::vector<size_t> neighs;
			for (size_t j = 0; j < temp.size(); ++j) {
				auto nei = getMO_HLC().getSubspaceInfo(temp[j]).m_sub_neighbors;
				for (auto jj : nei) {
					if (neighs.empty()) {
						neighs.push_back(jj);
					}
					else if (std::find(neighs.begin(),neighs.end(),jj)==neighs.end()&& std::find(temp.begin(), temp.end(), jj) == temp.end()) {
						neighs.push_back(jj);
					}
				}
			}
			temp_cluster.emplace_back(neighs);
			getMO_HLC().getClusters().emplace_back(temp_cluster);
		}
	}

	std::vector<std::vector<size_t>> SPMOEA1_8::clusterFrontSpace(const std::vector<size_t>& frontspace) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//标记
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space = i;
					head_cluster.push_back(frontspace[begin_space]);//先加入一个子空间
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
			while (count < select_flag.size()) {
				std::vector<size_t> temp;
				for (size_t j = 0; j < temp_cluster.size(); ++j) {
					size_t inx = temp_cluster[j];
					std::list<size_t> neighbors;
					getMO_HLC().subspaceTree().findNeighbor(inx, neighbors);
					for (size_t k = 0; k < frontspace.size(); ++k) {
						if (select_flag[k] == 0) {
							if (std::find(neighbors.begin(), neighbors.end(), frontspace[k]) != neighbors.end()) {
								////邻域子空间的连续性判断
								//if (subspaceLink(inx,frontspace[k])) {
								//	head_cluster.push_back(frontspace[k]);
								//	temp.push_back(frontspace[k]);
								//	select_flag[k] = 1;
								//	count++;
								//}
								head_cluster.push_back(frontspace[k]);
								temp.push_back(frontspace[k]);
								select_flag[k] = 1;
								count++;
							}
						}
					}
				}
				if (temp.empty()) {//没有新的邻域加入
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
		getMO_HLC().setFrontClusters(clustered);
		return clustered;
	}

	bool SPMOEA1_8::subspaceLink(size_t inx1, size_t inx2) {
		bool flag = false;
		//两个子空间中的前沿解的分布
		auto ind1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
		auto ind2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
		////到子空间边界的平均距离
		auto bound1 = getMO_HLC().subspaceTree().getBox(inx1);
		auto bound2 = getMO_HLC().subspaceTree().getBox(inx2);
		//std::map<size_t,Real> dim_dist1;
		//for (size_t i = 0; i < bound1.size(); ++i) {
		//	Real temp_sum = 0.;
		//	for (size_t j = 0; j < ind1.size(); ++j) {
		//		temp_sum += ind1[j]->variable()[i];
		//	}
		//	temp_sum /= ind1.size();
		//	Real min_dist = std::fabs(temp_sum - bound1[i].first) < std::fabs(temp_sum - bound1[i].second) ? std::fabs(temp_sum - bound1[i].first) : std::fabs(temp_sum - bound1[i].second);
		//	Real min_dist_ratio = min_dist / std::fabs(bound1[i].second - bound1[i].first);
		//	dim_dist1.insert(std::make_pair(i,min_dist_ratio));
		//}

		//std::map<size_t, Real> dim_dist2;
		//for (size_t i = 0; i < bound2.size(); ++i) {
		//	Real temp_sum = 0.;
		//	for (size_t j = 0; j < ind2.size(); ++j) {
		//		temp_sum += ind2[j]->variable()[i];
		//	}
		//	temp_sum /= ind2.size();
		//	Real min_dist = std::fabs(temp_sum - bound2[i].first) < std::fabs(temp_sum - bound2[i].second) ? std::fabs(temp_sum - bound2[i].first) : std::fabs(temp_sum - bound2[i].second);
		//	Real min_dist_ratio = min_dist / std::fabs(bound2[i].second - bound2[i].first);
		//	dim_dist2.insert(std::make_pair(i, min_dist_ratio));
		//}
		//for (size_t i = 0; i < dim_dist1.size();++i) {

		//}
		//计算两类点中的最短距离值
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < ind1.size(); ++i) {
			for (size_t j = 0; j < ind2.size(); ++j) {
				auto dist = euclideanDistance(ind1[i]->variable().vect().begin(), ind1[i]->variable().vect().end(), ind2[j]->variable().vect().begin());
				min_v = min_v < dist ? min_v : dist;
			}
		}
		//两类点中各自的点间最大值
		Real max_dist1 = 0.;
		Real max_dist2 = 0;
		std::vector<Real> dim_span1;
		std::vector<Real> dim_span2;
		//寻找每一维的切面
		for (size_t i = 0; i < bound1.size(); ++i) {
			Real min_v = (Real)INT16_MAX;
			Real max_v = -1. * INT16_MAX;
			for (size_t j = 0; j < ind1.size(); ++j) {
				if (min_v > ind1[j]->variable().vect()[i]) {
					min_v = ind1[j]->variable().vect()[i];
				}
				if (max_v < ind1[j]->variable().vect()[i]) {
					max_v = ind1[j]->variable().vect()[i];
				}
			}
			dim_span1.push_back(max_v - min_v);
		}
		for (size_t i = 0; i < dim_span1.size(); ++i) {
			max_dist1 += std::pow(dim_span1[i], 2);
		}
		for (size_t i = 0; i < bound2.size(); ++i) {
			Real min_v = (Real)INT16_MAX;
			Real max_v = -1. * INT16_MAX;
			for (size_t j = 0; j < ind2.size(); ++j) {
				if (min_v > ind2[j]->variable().vect()[i]) {
					min_v = ind2[j]->variable().vect()[i];
				}
				if (max_v < ind2[j]->variable().vect()[i]) {
					max_v = ind2[j]->variable().vect()[i];
				}
			}
			dim_span2.push_back(max_v - min_v);
		}
		for (size_t i = 0; i < dim_span2.size(); ++i) {
			max_dist2 += std::pow(dim_span2[i], 2);
		}
		max_dist1 = std::sqrt(max_dist1);
		max_dist2 = std::sqrt(max_dist2);
		if (min_v < max_dist1 / 3 || min_v < max_dist2 / 3) {
			flag = true;
		}

		return flag;
	}

	void SPMOEA1_8::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_8::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem *pro, Random *rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA1_8::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_8::findSplitDim(int inx, Problem *pro) {
		//根据在哪一维上平分两边的子空间前沿解的差额最大来确定划分的维度
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> front_ind_ratio;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			Real middle_v = (space_bound[j].first + space_bound[j].second) / 2;
			size_t upper_count = 0, lower_count = 0;
			for (size_t k = 0; k < front_sols.size(); ++k) {
				if (front_sols[k]->variable()[j] >= middle_v) {
					upper_count++;
				}
				else {
					lower_count++;
				}
			}
			size_t max_v = upper_count > lower_count ? upper_count : lower_count;
			size_t min_v = upper_count < lower_count ? upper_count : lower_count;
			if (max_v == 0) {
				max_v += 1;
			}
			front_ind_ratio.push_back((Real)min_v / max_v);
		}
		int dim = std::distance(front_ind_ratio.begin(), std::min_element(front_ind_ratio.begin(), front_ind_ratio.end()));

		//使用跨度确定分割的维度
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		return dim;
	}

	void SPMOEA1_8::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i) {
			objs.emplace_back(&pop[i]->objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i) {
			pop[i]->setFitness(rank[i]);
		}
	}

	void SPMOEA1_8::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		//使用archive计算性能指标
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			temp_pop.append(*m_archive[i]);
		}
		Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);
		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

	void SPMOEA1_8::initiObjSpace(Problem *pro) {

	}
}