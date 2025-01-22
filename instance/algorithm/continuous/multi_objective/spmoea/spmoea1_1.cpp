#include "spmoea1_1.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_1::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_1::run_() {
		initPop(m_problem.get(), this, m_random.get());
		while (!terminating()) {
			evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void SPMOEA1_1::record() {
		//std::vector<Real> entry;
		//entry.push_back(m_evaluations);

		//Population<Solution<>> temp_pop;
		//for (size_t i = 0; i < m_archive.size(); ++i) {
		//	temp_pop.append(*m_archive[i]);
		//}
		////IGD
		//Real IGD = m_problem->optimaBase()->invertGenDist(temp_pop);
		//entry.push_back(IGD);
		////HV
		////先得到PF采样点中的最大值作为ref_point
		//auto& refPoint = getMaxRefPoint();
		//Real HV = hypervolumePop(temp_pop, refPoint, m_problem.get(), m_random.get());
		//entry.push_back(HV);
		//dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_1::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(2 * getPop().size() + 1);//第二个为子代
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j]);
				}
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i + 1].push_back(&getPop()[i].getOffspring()[j].phenotype());
				}
			}
			auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(his_sols[i].get());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA1_1::initiVarSpace(Problem* pro) {
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

	void SPMOEA1_1::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//在各个子空间生成一个子种群
		size_t number_objectives = CAST_CONOP(pro)->numberObjectives();
		size_t num_vars = CAST_CONOP(pro)->numberVariables();
		if (number_objectives == 2) {
			if (num_vars == 2) {
				setArchiveNum(20);
			}
			else if (num_vars == 5) {
				setArchiveNum(50);
			}
			else if (num_vars == 10) {
				setArchiveNum(100);
			}
			else {
				setArchiveNum(100);
			}
		}
		else if (number_objectives == 3) {
			if (num_vars == 3) {
				setArchiveNum(100);
			}
			else if (num_vars == 5) {
				setArchiveNum(100);
			}
			else if (num_vars == 10) {
				setArchiveNum(100);
			}
			else {
				setArchiveNum(200);
			}
		}
		else {
			setArchiveNum(500);
		}
		size_t pop_num = getPopsize();
		Population<Solution<>> new_pop;
		//全局种群
		SPMOEA_pop temp_pop(pop_num, pro);
		temp_pop.initialize(pro, rnd);
		std::vector<std::vector<Real>> sols;
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t j = 0; j < pop_num; ++j) {
			temp_pop[j].evaluate(pro, alg);
			sols.emplace_back(temp_pop[j].variable().vect());
		}
		SPMOEA::NDSort(temp_pop);
		//初始化子代
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			new_pop.append(temp_pop[j]);
		}
		temp_pop.setPopState("active");
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

		//更新子空间信息
		for (size_t i = 0; i < getPop().size(); ++i) {
			Population<Solution<>> temp_off_pop;
			for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); ++j) {
				temp_off_pop.append(getPop()[i].getOffspring()[j]);
			}
			SPMOEA::updateSubspaceFrontSol(temp_off_pop, pro, rnd);
		}
		updateFrontSpace();
		updateEE(pop_num, 0);

		updateFrontRegionLinkSpace(pro,rnd);
		//根据子空间划分聚类
		clusterSubspace();
		setAddNumSpace(0);
		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < getFrontSpace().size(); ++i) {
			auto v = getMO_HLC().subspaceTree().getBoxVolume(getFrontSpace()[i]);
			front_ratio += (v / total_volume);
		}
		setFrontSpaceRatio(front_ratio);
		m_divide_granularity = 20 - 10. / 8 * ((num_vars >= 10 ? 10 : num_vars) - 2);
		//加入历史演化轨迹
		for (size_t j = 0; j < getPop()[0].size(); ++j) {
			std::vector<std::shared_ptr<Solution<>>> temp_ind;
			temp_ind.emplace_back(std::make_shared<Solution<>>(getPop()[0][j]));
			getHisEvolveLocus().emplace_back(temp_ind);
		}
		//算子比例
		std::vector<Real> operator_ratio(2, 1.);
		operator_ratio[1] = 0.;
		getOperatorRatio().emplace_back(operator_ratio);
		//个体所在的子空间
		std::vector<size_t> ind_space;
		for (size_t j = 0; j < getPop()[0].size(); ++j) {
			auto sol = getPop()[0][j].variable().vect();
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			if (std::find(ind_space.begin(), ind_space.end(), space_inx) == ind_space.end()) {
				ind_space.push_back(space_inx);
			}
		}
		std::vector<std::vector<size_t>> plot_spaces;
		auto front_clusters = getFrontRegionLinkSpace();
		auto front_spaces = getFrontSpace();
		for (size_t j = 0; j < front_clusters.size(); ++j) {
			std::vector<size_t> temp;
			for (size_t k = 0; k < ind_space.size(); ++k) {
				if (std::find(front_clusters[j].begin(), front_clusters[j].end(), ind_space[k]) != front_clusters[j].end()) {
					temp.push_back(ind_space[k]);
				}
			}
			plot_spaces.emplace_back(temp);
			//plot_spaces.emplace_back(front_clusters[j]);
		}
		for (size_t k = 0; k < ind_space.size(); ++k) {
			if (std::find(front_spaces.begin(), front_spaces.end(), ind_space[k]) == front_spaces.end()) {
				std::vector<size_t> temp;
				temp.push_back(ind_space[k]);
				plot_spaces.emplace_back(temp);
			}
		}
		getIndSpaces() = plot_spaces;

		//SPMOEA::recordMetrics(pro, alg);
		SPMOEA::record();
	}

	int SPMOEA1_1::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		/********************************************************************************
							   根据子种群前沿子空间的体积分配搜索资源
		********************************************************************************/
		//根据算子类型和子连通区域长度分配计算资源
		getInteractiveSols().clear();
		std::vector<size_t> assign_pop_resource;
		size_t switch_period = 3;
		PopResourceAssign(assign_pop_resource, switch_period, pro);
		/********************************************************************************
											 子种群演化
		********************************************************************************/
		//1.子代中探索与开发的比例；2.子代中探索与开发的方式
		//根据子空间前沿解的线性化程度和体积大小决定子代生成方式
		//若非线性强，
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			std::vector<int> interactive_type;//可以设置不同连通子空间的交互类型
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (m_divide_iteration % switch_period == 0) {
					interactive_type.push_back(1);
				}
				else {
					interactive_type.push_back(2);
				}
			}
			generateOffspring(pro, alg, rnd, assign_pop_resource, interactive_type);
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
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}
		SPMOEA::NDSort(offspring_pop);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		updateSubspaceFrontSol(offspring_pop, pro, rnd);//子代更新子空间前沿解和代表解
		auto pre_front_space = getFrontSpace();
		updateFrontSpace();//根据前沿解更新子空间信息
		/**********************************************************************************
					       检测是否细分子空间:子空间前沿解的线性化程度
		***********************************************************************************/
		auto pre_num_spaces = getMO_HLC().numSubspace();
		size_t divide_fre = 10;
		if (m_divide_iteration % divide_fre == 0) {
			//检索前沿子空间是否可分
			auto cur_front_spaces = getFrontSpace();
			for (size_t i = 0; i < cur_front_spaces.size(); ++i) {
				//子空间解的线性判别和子连通空间的线性判别的阈值设置
				int flag = subspaceSeperable(cur_front_spaces[i], pro);
				if (flag) {
					//子空间划分维度的选择还需要优化，最好能够将非线性的分布划分为两个线性的分布
					spaceSubdivision(cur_front_spaces[i],pro, rnd);
				}
				else {
					getMO_HLC().getSubspaceInfo(cur_front_spaces[i]).m_linear_flag = 1;
				}
			}
			updateFrontSpace();
			//更新子空间邻域
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
			}
		}
		auto cur_num_spaces = getMO_HLC().numSubspace();
		setAddNumSpace(cur_num_spaces - pre_num_spaces);
		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		updateFrontRegionLinkSpace(pro,rnd);

		//测试分散子空间是否满足线性关系
		auto &clusters = getMO_HLC().getClusters();
		std::vector<size_t> single_space;
		for (size_t i = 0; i < clusters.size(); ++i) {
			if (clusters[i][0].size() == 1) {
				single_space.push_back(clusters[i][0][0]);
			}
		}
		auto front_sp = getFrontSubspace();
		auto temp_cluster = clusterRegionFrontSpace(single_space);
		for (size_t i = 0; i < temp_cluster.size(); ++i) {
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < temp_cluster[i].size(); ++j) {
				auto& front_sols1 = getMO_HLC().getSubspaceInfo(temp_cluster[i][j]).m_front_sol_in_subspace;
				for (size_t k = 0; k < front_sols1.size(); ++k) {
					all_vars.emplace_back(front_sols1[k]->variable().vect());
				}
			}
			if (all_vars.size() > CAST_CONOP(pro)->numberVariables()) {
				std::vector<std::vector<Real>> temp_vars = all_vars;
				auto bound = CAST_CONOP(pro)->boundary();
				//bound取点的上下界
				for (size_t i = 0; i < bound.size(); ++i) {
					bound[i].first = CAST_CONOP(pro)->boundary()[i].second;
					bound[i].second = CAST_CONOP(pro)->boundary()[i].first;
					for (size_t j = 0; j < temp_vars.size(); ++j) {
						if (bound[i].first > temp_vars[j][i]) {
							bound[i].first = temp_vars[j][i];
						}
						if (bound[i].second < temp_vars[j][i]) {
							bound[i].second = temp_vars[j][i];
						}
					}
				}
				dataNormalizeInBound(temp_vars, bound);
				//归一化数据后，求这组点的线性组合，采用奇异值分解求系数矩阵
				Eigen::MatrixXd Apoints = Eigen::MatrixXd::Random(temp_vars.size(), bound.size());
				Eigen::MatrixXd bconstant = Eigen::MatrixXd::Ones(temp_vars.size(), 1);
				Eigen::MatrixXd Xk1, Xk2;
				for (size_t i = 0; i < temp_vars.size(); ++i) {
					for (size_t j = 0; j < bound.size(); ++j) {
						Apoints(i, j) = temp_vars[i][j];
					}
				}
				Xk1 = Apoints.colPivHouseholderQr().solve(bconstant);
				//用得到的线性组合反算误差
				std::vector<Real> cal_error(temp_vars.size(), 0.);
				Real mean_error = 0.;
				for (size_t i = 0; i < temp_vars.size(); ++i) {
					auto temp = Apoints.row(i) * Xk1;
					Real error = std::fabs(temp.value() - 1);
					cal_error[i] = error;
					mean_error += error;
				}
				mean_error /= cal_error.size();
				size_t count = 0;
				for (size_t i = 0; i < cal_error.size(); ++i) {
					if (cal_error[i] > mean_error) {//误差百分比
						count++;
					}
				}
				Real error_ratio = (Real)count / cal_error.size();
				std::cout << "the error ratio is " << error_ratio << std::endl;
			}
		}
		
		//更新线性邻域子空间
		auto front_link_spaces = getFrontRegionLinkSpace();
		for (size_t i = 0;i< front_link_spaces.size(); ++i) {
			if (front_link_spaces[i].size() > 1) {
				for (size_t j = 0; j < front_link_spaces[i].size(); ++j) {
					for (size_t k = 0; k < front_link_spaces[i].size(); ++k) {
						if (j != k) {
							if (ifLinearLinkSpace(front_link_spaces[i][j], front_link_spaces[i][k])) {
								getMO_HLC().getSubspaceInfo(front_link_spaces[i][j]).m_linear_neigh_space_inx.push_back(front_link_spaces[i][k]);
							}
						}
						else {
							getMO_HLC().getSubspaceInfo(front_link_spaces[i][j]).m_linear_neigh_space_inx.push_back(front_link_spaces[i][j]);
						}
					}
				}
			}
			else {
				getMO_HLC().getSubspaceInfo(front_link_spaces[i][0]).m_linear_neigh_space_inx.push_back(front_link_spaces[i][0]);
			}
		}
		setFrontLastGens(1);
		clusterSubspace();
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		//SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < getFrontSpace().size(); ++i) {
			auto v = getMO_HLC().subspaceTree().getBoxVolume(getFrontSpace()[i]);
			front_ratio += (v / total_volume);
			front_count++;
		}
		setFrontSpaceRatio(front_ratio);
		/********************************************************************************************
										 种群更新，子种群内部进行淘汰选择
		*********************************************************************************************/
		updateInteractiveSols();
		//采用局部比较的方式，选择个体
		//在进行目标空间排序选择时，考虑解在决策空间的距离
		//multiObjSelection(pro);
		//sparseSelection(pro);
		//localSelection(pro);
		//localCrowdSelection(pro,rnd);
		//localCrowdSelection2(pro, rnd);
		ensembleSelection(getPop()[0].size(), pro, rnd);

		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		//SPMOEA::recordMetrics(pro, alg);
		SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}


	void SPMOEA1_1::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA1_1::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_1::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
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

	void SPMOEA1_1::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		auto front_link_spaces = getFrontRegionLinkSpace();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		size_t total_size = 0;
		//根据连通子空间的前沿个体比例分配个体数
		std::vector<size_t> front_num_in_clusters;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			size_t count = 0;
			for (size_t j = 0; j < front_link_spaces[i].size(); ++j) {
				count += getMO_HLC().getSubspaceInfo(front_link_spaces[i][j]).m_front_sol_in_subspace.size();
			}
			total_size += count;
			front_num_in_clusters.push_back(count);
		}
		size_t count = 0;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			size_t temp_size = (size_t)std::ceil(max_pop_size*front_num_in_clusters[i] / (Real)total_size);
			assign_pop_resource.push_back(temp_size);
			count += temp_size;
		}
		while (count > max_pop_size) {
			auto inx=std::distance(assign_pop_resource.begin(), std::max_element(assign_pop_resource.begin(), assign_pop_resource.end()));
			assign_pop_resource[inx]--;
			count--;
		}
		
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_1::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
		size_t num_exploit = 0;
		auto search_bound = CAST_CONOP(pro)->boundary();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		//根据子连通集采样
		auto front_clusters = getFrontRegionLinkSpace();
		auto front_spaces = getFrontSpace();
		for (size_t i = 0; i < front_clusters.size(); ++i) {
			//连通子空间内的个体形成种群
			size_t kk = 2;
			std::vector<std::vector<Real>> all_off;
			SPMOEA_pop temp_pop(0, pro);//临时的交互种群
			for (size_t j = 0; j < getPop()[0].size(); ++j) {
				auto& sol = getPop()[0][j].variable().vect();
				auto space = getMO_HLC().subspaceTree().getRegionIdx(sol);
				if (std::find(front_clusters[i].begin(), front_clusters[i].end(), space) != front_clusters[i].end()) {
					temp_pop.append(getPop()[0][j]);
				}
			}
			if (temp_pop.size() < 5) {
				//不够的在子连通集内变异
				std::vector<std::shared_ptr<Solution<>>> temp_his_sols;
				for (size_t j = 0; j < front_clusters[i].size(); ++j) {
					auto& his_sols = getMO_HLC().getSubspaceInfo(front_clusters[i][j]).m_history_inds;
					for (size_t k = 0; k < his_sols.size(); ++k) {
						temp_his_sols.emplace_back(his_sols[k]);
					}
				}
				if (temp_his_sols.size() < 5) {
					std::vector<std::vector<Real>> offs;
					for (size_t j = 0; j < pop_resource[i]; ++j) {
						size_t inx= (size_t)std::floor(front_clusters[i].size() * rnd->uniform.next());
						auto& front_sol = getMO_HLC().getSubspaceInfo(front_clusters[i][inx]).m_subspace_front_sol;
						size_t idx= (size_t)std::floor(front_sol.size() * rnd->uniform.next());
						offs=sampleInRange(*front_sol[idx], 1, pro, alg, rnd);
						for (auto off : offs) {
							all_off.emplace_back(off);
						}
					}
				}
				else {
					size_t count = 0;
					while (temp_pop.size() < 5) {
						size_t inx = (size_t)std::floor(temp_his_sols.size() * rnd->uniform.next());
						auto& sol = temp_his_sols[inx]->variable().vect();
						bool flag = true;
						for (size_t j = 0; j < temp_pop.size(); ++j) {
							if (ifSame(sol, temp_pop[j].variable().vect())) {
								flag = false;
								break;
							}
						}
						if (flag) {
							temp_pop.append(*temp_his_sols[inx]);
						}
						count++;
						if (count > 10000) {
							break;
						}
					}
					if (temp_pop.size() < 5) {
						std::vector<std::vector<Real>> offs;
						for (size_t j = 0; j < pop_resource[i]; ++j) {
							size_t inx = (size_t)std::floor(front_clusters[i].size() * rnd->uniform.next());
							auto& front_sol = getMO_HLC().getSubspaceInfo(front_clusters[i][inx]).m_subspace_front_sol;
							size_t idx = (size_t)std::floor(front_sol.size() * rnd->uniform.next());
							offs = sampleInRange(*front_sol[idx], 1, pro, alg, rnd);
							for (auto off : offs) {
								all_off.emplace_back(off);
							}
						}
					}
					else {
						all_off = sampleByDE(temp_pop, search_bound, pop_resource[i], kk, pro, alg, rnd);
					}
					
				}
			}
			else {
				//种群在连通子空间交互
				all_off = sampleByDE(temp_pop, search_bound, pop_resource[i], kk, pro, alg, rnd);
			}
		}
		auto& interactive_sols = getInteractiveSols();
		for (size_t k = 0; k < interactive_sols.size(); ++k) {
			auto ind = *interactive_sols[k].back();
			getPop()[0].getOffspring()[k].variable() = ind.variable();
			getPop()[0].getOffspring()[k].objective() = ind.objective();
			getPop()[0].getOffspring()[k].setCounter(0);
		}
		
		updateEE(0, num_exploit);
	}

	void SPMOEA1_1::updateFrontRegionLinkSpace(Problem* pro, Random* rnd) {
		//先找到前沿子空间
		std::vector<size_t> front_space = getFrontSpace();
		//auto link_space = clusterRegionFrontSpace(front_space);
		auto link_space = linearClusterFrontSpace(front_space,pro,rnd);//线性子空间聚类
		getFrontRegionLinkSpace() = link_space;
	}

	bool SPMOEA1_1::spaceSubdivision(size_t space_inx,Problem* pro, Random* rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		//细分子空间
		auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(space_inx);
		size_t min_num = m_divide_granularity;
		if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
			//先得到子空间历史解用于更新细分的子空间的信息
			/*int dim = findSplitDim(space_inx, pro);
			auto& space_bound = getMO_HLC().subspaceTree().getBox(space_inx);
			Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
			splitSpace(space_inx, 2, dim, pos, false, pro, rnd);*/
			auto divide_info = findSplitDimAndPos(space_inx,pro);
			splitSpace(space_inx, 2, divide_info.first, divide_info.second, false, pro, rnd);
		}
		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}
		return flag;
	}

	bool SPMOEA1_1::spaceSubdivision(Problem* pro, Random* rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		auto front_spaces = getFrontSpace();
		//细分子空间
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			size_t min_num = m_divide_granularity;
			//size_t min_num = getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity;
			// min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
			//min_num = 20;
			if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
				//先得到子空间历史解用于更新细分的子空间的信息
				int dim = findSplitDim(sp, pro);
				auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
				Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
				splitSpace(sp, 2, dim, pos, false, pro, rnd);
			}
		}

		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}

		//子空间体积占比
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		return flag;
	}

	void SPMOEA1_1::clusterSubspace() {
		//每次聚类只加入一层，允许类间有重叠
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
		}
		getMO_HLC().getClusters().clear();
		/*for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					auto inx = getMO_HLC().getClusters()[i][j][k];
					getMO_HLC().getSubspaceInfo(inx).idx_cluster.clear();
				}
			}
		}*/
		auto& front_clusters = getFrontRegionLinkSpace();
		
		for (size_t i = 0; i < front_clusters.size(); ++i) {
			std::vector<std::vector<size_t>> clusters;
			clusters.emplace_back(front_clusters[i]);
			getMO_HLC().getClusters().emplace_back(clusters);
		}
		for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					auto inx = getMO_HLC().getClusters()[i][j][k];
					getMO_HLC().getSubspaceInfo(inx).idx_cluster.push_back(i);
				}
			}
		}
	}

	std::vector<std::vector<size_t>> SPMOEA1_1::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_1::subspaceLink(size_t inx1, size_t inx2) {
		bool flag = false;
		//两个子空间中的前沿解的分布
		auto ind1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
		auto ind2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
		////到子空间边界的平均距离
		auto bound1 = getMO_HLC().subspaceTree().getBox(inx1);
		auto bound2 = getMO_HLC().subspaceTree().getBox(inx2);
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

	void SPMOEA1_1::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	//return 1: seperable
	int SPMOEA1_1::subspaceSeperable(size_t inx, Problem* pro) {
		// 线性判定和超平面判定
        // 线性判定采用角度偏差，超平面判定采用线性方程组误差
		//obtain front sols in inx subspace
		std::vector<std::vector<Real>> all_vars;
		size_t var_dim = CAST_CONOP(pro)->numberVariables();
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_front_sol_in_subspace;
		for (size_t i = 0; i < front_sols.size(); ++i) {
			all_vars.emplace_back(front_sols[i]->variable().vect());
		}
		if (all_vars.size() <= var_dim) {//点的个数小于维数不可分
			return 0;
		}
		else {
			std::vector<std::vector<Real>> temp_vars = all_vars;
			auto bound = getMO_HLC().subspaceTree().getBox(inx);
			//bound取点的上下界
			for (size_t i = 0; i < bound.size(); ++i) {
				bound[i].first = CAST_CONOP(pro)->boundary()[i].second;
				bound[i].second = CAST_CONOP(pro)->boundary()[i].first;
				for (size_t j = 0; j < temp_vars.size(); ++j) {
					if (bound[i].first > temp_vars[j][i]) {
						bound[i].first = temp_vars[j][i];
					}
					if (bound[i].second < temp_vars[j][i]) {
						bound[i].second = temp_vars[j][i];
					}
				}
			}
			dataNormalizeInBound(temp_vars, bound);
			////归一化数据后，求这组点的线性组合，采用奇异值分解求系数矩阵
			//Eigen::MatrixXd Apoints = Eigen::MatrixXd::Random(temp_vars.size(), var_dim);
			//Eigen::MatrixXd bconstant = Eigen::MatrixXd::Ones(temp_vars.size(), 1);
			//Eigen::MatrixXd Xk1,Xk2;
			//for (size_t i = 0; i < temp_vars.size(); ++i) {
			//	for (size_t j = 0; j < var_dim; ++j) {
			//		Apoints(i, j) = temp_vars[i][j];
			//	}
			//}
			//Xk1 = Apoints.colPivHouseholderQr().solve(bconstant);
			////Xk2 = Apoints.jacobiSvd(Eigen::ComputeThinV | Eigen::ComputeThinU).solve(bconstant);
			////用得到的线性组合反算误差
			//std::vector<Real> cal_error(temp_vars.size(), 0.);
			////std::vector<Real> cal_error2(temp_vars.size(), 0.);
			//Real mean_error = 0.;
			//for (size_t i = 0; i < temp_vars.size(); ++i) {
			//	auto temp = Apoints.row(i) * Xk1;
			//	Real error = std::fabs(temp.value() - 1);
			//	cal_error[i] = error;
			//	mean_error += error;
			//	/*auto temp2 = Apoints.row(i) * Xk2;
			//	Real error2 = std::fabs(temp2.value() - 1);
			//	cal_error2[i] = error2;*/
			//}
			//mean_error/=cal_error.size();
			
			//统计向量夹角
			Real max_dist = 0;
			std::pair<size_t, size_t> max_inx;
			for (size_t i = 0; i < temp_vars.size(); ++i) {
				for (size_t j = i + 1; j < temp_vars.size(); ++j) {
					Real temp_dist = euclideanDistance(temp_vars[i].begin(),temp_vars[i].end(),temp_vars[j].begin());
					if (temp_dist > max_dist) {
						max_dist=temp_dist;
						max_inx.first = i;
						max_inx.second = j;
					} 
				}
			}
			//计算角点位置
			auto corner = temp_vars[max_inx.first];
			/*for (size_t i = 0; i < var_dim; ++i) {
				corner[i] = 2 * corner[i] - temp_vars[max_inx.second][i];
			}*/
			//统计所有点与向量之间的距离
			std::vector<Real> angles;
			std::vector<Real> vertical_dist;
			for (size_t i = 0; i < temp_vars.size(); ++i) {
				Real ang = vectorAngle(corner, temp_vars[max_inx.second],temp_vars[i]);
				Real l = euclideanDistance(temp_vars[i].begin(),temp_vars[i].end(),corner.begin());
				vertical_dist.push_back(l*std::sin(ang));
				angles.push_back(ang/OFEC_PI*180);
			}
			size_t count = 0;
			for (size_t i = 0; i < vertical_dist.size(); ++i) {
				if (vertical_dist[i] > max_dist*0.5) {//误差百分比
					count++;
				}
			}
			if (count > angles.size() / 2) {
				return 1;
			}
			//for (size_t i = 0; i < cal_error.size(); ++i) {
			//	if (cal_error[i] > mean_error) {//误差百分比
			//		count++;
			//	}
			//}
			//if (count > cal_error.size() / 3) {
			//	return 1;
			//}
			
			////奇异值分解，根据特征值判断线性程度
			//std::vector<Real> vector_data;
			//for (size_t i = 0; i < temp_vars.size(); ++i) {
			//	for (size_t j = 0; j < var_dim; ++j) {
			//		vector_data.push_back(temp_vars[i][j]);
			//	}
			//}
			//Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix_data(vector_data.data(),temp_vars.size(), var_dim);
			//Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix_data, Eigen::ComputeThinV | Eigen::ComputeThinU);
			//Eigen::MatrixXd singular_values = svd.singularValues();
			////比较奇异值的大小，确定线性程度和维度
			//std::vector<Real> qr_value;
			//for (size_t i = 0; i < singular_values.size(); ++i) {
			//	qr_value.push_back(singular_values(i));
			//}
			////通过最后一个奇异值的大小与前一个奇异值的大小判断
			////或者有效的特征值的数量小于维数
			//for (size_t i = qr_value.size()-1; i>1; --i) {
			//	Real qr1 = qr_value[i];
			//	Real qr2 = qr_value[i-1];
			//	if (qr1 / qr2 < 0.1||qr1<10e-5) {
			//		return 1;
			//	}
			//}
		}
		return 0;
	}

	void SPMOEA1_1::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA1_1::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_1::findSplitDim(int inx, Problem* pro) {
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

	

	void SPMOEA1_1::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_1::recordMetrics(Problem* pro, Algorithm* alg) {
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

	void SPMOEA1_1::initiObjSpace(Problem* pro) {

	}
}