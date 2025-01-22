#include "spmoea2_1.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA2_1::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA2_1::run_() {
		initPop(m_problem.get(), this, m_random.get());
		//#ifdef OFEC_DEMO
		//		updateBuffer();
		//#endif
		while (!terminating()) {
			evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void SPMOEA2_1::record() {
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
	void SPMOEA2_1::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(2 * getPop().size() + 1);//第二个为子代
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j].phenotype());
				}
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i + 1].push_back(&getPop()[i].getOffspring()[j].phenotype());
				}
			}
			auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(&his_sols[i]->phenotype());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA2_1::initiVarSpace(Problem* pro) {
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

	void SPMOEA2_1::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
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
		//size_t space_num = getMO_HLC().numSubspace();//初始化子空间数量
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
		//temp_pop.setRate(getCr(), getMr());
		//temp_pop.setEta(getCeta(), getMeta());
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

		updateFrontRegionLinkSpace(pro, rnd);
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
		m_divide_granularity = getSplitGranularity();
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

	int SPMOEA2_1::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
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
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			bool m_assign_resource = true;//子种群是否进行资源分配
			if (!m_assign_resource) {
				for (auto& i : assign_pop_resource) {
					i = getPopsize();
				}
			}
			std::vector<int> interactive_type;
			auto front_clusters = getFrontRegionLinkSpace();
			auto front_spaces = getFrontSpace();
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (m_divide_iteration % switch_period == 0) {
					interactive_type.push_back(3);
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
				/*tag = getPop()[i].getOffspring()[j].evaluate(pro, alg);
				if (tag != EvaluationTag::kNormalEval)
					break;*/
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}
		SPMOEA::NDSort(offspring_pop);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		updateSubspaceFrontSol(offspring_pop, pro, rnd);//子代更新子空间前沿解和代表解
		auto pre_front_space = getFrontSpace();
		updateFrontSpace();//使用子种群子代更新子空间信息
		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		auto pre_num_spaces = getMO_HLC().numSubspace();
		if (m_divide_iteration % 1 == 0) {
			auto divide_flag = spaceSubdivision(pro, rnd);
			if (!divide_flag) {
				m_stage_last_time++;
				if (m_stage_last_time > 2 && m_divide_granularity < 50) {
					m_divide_granularity++;
					m_stage_last_time = 0;
				}
				//spaceSubdivision(pro, rnd);
			}
			else {
				m_stage_last_time = 0;
				//更新子空间邻域
				for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
					getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
				}
				updateFrontSpace();
			}
		}
		/*auto flag = ifFrontChanged(pre_front_space);
		if (!flag && getFrontLastGens() >= 3) {
			auto divide_flag=spaceSubdivision(pro, rnd);
			if (!divide_flag) {
				m_divide_granularity++;
				m_stage_converge = true;
				//spaceSubdivision(pro, rnd);
			}
			updateFrontSpace();
			setFrontLastGens(1);
			clusterSubspace();
		}*/
		auto cur_num_spaces = getMO_HLC().numSubspace();
		setAddNumSpace(cur_num_spaces - pre_num_spaces);
		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		updateFrontRegionLinkSpace(pro, rnd);
		setFrontLastGens(1);
		clusterSubspace();
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		//SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		////检测子空间是否包含不同的PS片段
		////auto front_sp = getFrontSpace();
		//auto front_sp = getFrontSpace();
		////真实PS所在的子空间
		//std::vector<size_t> real_front_subspaces;
		////具有多段的子空间数目
		//std::vector<size_t> multi_segment_subspaces;
		//for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//	auto& box = getMO_HLC().subspaceTree().getBox(i);
		//	//沿着第一维构造解
		//	size_t num = 500;
		//	std::vector<size_t> flag(num - 1, 1);
		//	Real delta = (box[0].second - box[0].first) / num;
		//	for (size_t j = 0; j < num - 1; ++j) {
		//		VariableVector<> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
		//		auto sol = CAST_CONOP(pro)->createSolution(s);
		//		//判断sol是否在box内
		//		for (size_t k = 0; k < box.size(); ++k) {
		//			if (dynamic_cast<const VariableVector<>&>(sol->variableBase())[k] > box[k].second || dynamic_cast<const VariableVector<>&>(sol->variableBase())[k] < box[k].first) {
		//				flag[j] = 0;
		//				break;
		//			}
		//		}
		//	}
		//	//统计段数
		//	size_t segment_time = 0;
		//	for (size_t j = 1; j < flag.size(); ++j) {
		//		if (flag[0] == 0) {
		//			if (flag[j - 1] == 0 && flag[j] == 1) {
		//				segment_time++;
		//			}
		//		}
		//		else {
		//			if (flag[j - 1] == 1 && flag[j] == 0) {
		//				segment_time++;
		//			}
		//		}
		//	}
		//	if (segment_time > 0) {
		//		real_front_subspaces.push_back(i);
		//		if (segment_time > 1) {
		//			multi_segment_subspaces.push_back(i);
		//		}
		//	}
		//	//std::cout << "the number of segment of the " << i << " th subspace is: " << segment_time << std::endl;
		//	//1、看子空间的非线性程度

		//	//2、看子空间内的解是否可分为多个类，采用dbscan聚类

		//}
		////真实子空间与逼近子空间的差异
		//std::vector<size_t> match_subspace;
		//std::vector<size_t> error_subspace;
		//std::vector<size_t> loss_subspace;
		//std::vector<size_t> match_multi_seg_subspace;
		//for (size_t i = 0; i < front_sp.size(); ++i) {
		//	if (std::find(real_front_subspaces.begin(), real_front_subspaces.end(), front_sp[i]) != real_front_subspaces.end()) {
		//		match_subspace.push_back(front_sp[i]);
		//	}
		//	else {
		//		error_subspace.push_back(front_sp[i]);
		//	}
		//}
		//for (size_t i = 0; i < real_front_subspaces.size(); ++i) {
		//	if (std::find(front_sp.begin(), front_sp.end(), real_front_subspaces[i]) == front_sp.end()) {
		//		loss_subspace.push_back(real_front_subspaces[i]);
		//	}
		//}
		//for (size_t i = 0; i < match_subspace.size(); ++i) {
		//	if (std::find(multi_segment_subspaces.begin(), multi_segment_subspaces.end(), match_subspace[i]) == multi_segment_subspaces.end()) {
		//		match_multi_seg_subspace.push_back(match_subspace[i]);
		//	}
		//}
		//std::cout << "the number of total   subspace is: " << getMO_HLC().numSubspace() << std::endl;
		//std::cout << "the number of real-PS subspace is: " << real_front_subspaces.size() << std::endl;
		//std::cout << "the number of mul-seg subspace is: " << multi_segment_subspaces.size() << std::endl;
		//std::cout << "the number of matched subspace is: " << match_subspace.size() << std::endl;
		//std::cout << "match multiseg subspace number is: " << match_multi_seg_subspace.size() << std::endl;
		//std::cout << "the number of errored subspace is: " << error_subspace.size() << std::endl;
		//std::cout << "the number of lossed  subspace is: " << loss_subspace.size() << std::endl;


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


	void SPMOEA2_1::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA2_1::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA2_1::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
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

	void SPMOEA2_1::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		//updateFrontRegionLinkSpace();
		if (m_stage_last_time > 0) {
			auto front_link_spaces = getFrontRegionLinkSpace();
			size_t M = CAST_CONOP(pro)->numberObjectives();
			size_t total_size = 0;
			for (size_t i = 0; i < front_link_spaces.size(); ++i) {
				total_size += (1 * front_link_spaces[i].size());
			}
			assign_pop_resource.push_back(max_pop_size);
		}
		else {
			assign_pop_resource.push_back(max_pop_size);
		}
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA2_1::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
		size_t num_exploit = 0;
		auto search_bound = CAST_CONOP(pro)->boundary();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		//前面的种群在自己的连通集之间扩展
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t exploit_num = pop_resource[i];
			num_exploit += exploit_num;
			//子代个体数量处理
			while (getPop()[i].getOffspring().size() - getPop()[i].size() > pop_resource[i]) {
				getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
			}
			while (getPop()[i].getOffspring().size() - getPop()[i].size() < pop_resource[i]) {
				Solution<> temp_ind(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
				getPop()[i].getOffspring().append(temp_ind);
			}
			auto bound = CAST_CONOP(pro)->boundary();
			auto front_clusters = getFrontRegionLinkSpace();
			auto front_spaces = getFrontSpace();
			size_t total_spaces = getMO_HLC().numSubspace();
			std::vector<size_t> space_fre;
			for (size_t j = 0; j < total_spaces; ++j) {
				space_fre.push_back(getMO_HLC().getSubspaceInfo(j).m_sub_freq);
			}
			std::vector<Real> operator_ratio(2, 0.);
			std::vector<std::vector<Real>> all_off;
			//在子区域前排个体间开发,开发方式选择与邻域个体之间交互
			std::vector<size_t> ind_order;
			if (type[i] == 1) {
				//auto all_off = sampleByGA(getPop()[i], exploit_num, pro, rnd);
				size_t kk = 2;
				all_off = sampleByDE(getPop()[i], bound, exploit_num, kk, pro, alg, rnd);
			}
			else if (type[i] == 2) {
				//基于每个个体采样
				size_t kk = 2;
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					auto& sol1 = getPop()[i][j].variable().vect();
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol1);
					auto& front_sols = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
					bool flag = false;
					for (size_t k = 0; k < front_sols.size(); ++k) {
						auto& sol2 = front_sols[k]->variable().vect();
						if (ifSame(sol1, sol2)) {
							flag = true;
							//邻域子空间内的个体交互
							SPMOEA_pop temp_pop(0, pro);
							//找出其邻域子空间
							std::vector<size_t> nei;
							nei.push_back(space_inx);
							auto neigh = getMO_HLC().getSubspaceInfo(space_inx).m_sub_neighbors;
							for (auto jj : neigh) {
								nei.push_back(jj);
							}
							//将当前种群在邻域的个体加入
							size_t temp_inx = 0;
							for (size_t p = 0; p < getPop()[i].size(); ++p) {
								auto& sol = getPop()[i][p].variable().vect();
								auto space = getMO_HLC().subspaceTree().getRegionIdx(sol);
								if (std::find(nei.begin(), nei.end(), space) != nei.end()) {
									temp_pop.append(getPop()[i][p]);
								}
								if (ifSame(sol1, sol)) {
									temp_inx = temp_pop.size() - 1;
								}
							}
							if (temp_pop.size() < 5) {
								//看邻域历史解是否够
								size_t his_num = 0;
								for (size_t p = 0; p < nei.size(); ++p) {
									his_num += getMO_HLC().getSubspaceInfo(nei[p]).m_history_inds.size();
								}
								if (his_num >= 5) {
									while (temp_pop.size() < 5) {
										for (size_t p = 0; p < nei.size(); ++p) {
											auto& his_sol = getMO_HLC().getSubspaceInfo(nei[p]).m_history_inds;
											for (size_t q = 0; q < his_sol.size(); ++q) {
												auto& temp_sol = his_sol[q]->variable().vect();
												if (!ifSame(sol1, temp_sol)) {
													temp_pop.append(*his_sol[q]);
												}
												if (temp_pop.size() >= 5) {
													break;
												}
											}
											if (temp_pop.size() >= 5) {
												break;
											}
										}
									}
								}
								else {
									//加入非邻域个体
									for (size_t p = 0; p < getPop()[i].size(); ++p) {
										auto& sol = getPop()[i][p].variable().vect();
										auto space = getMO_HLC().subspaceTree().getRegionIdx(sol);
										if (std::find(nei.begin(), nei.end(), space) == nei.end()) {
											temp_pop.append(getPop()[i][p]);
										}
										if (temp_pop.size() >= 5) {
											break;
										}
									}
								}
							}
							auto off = sampleByDE(temp_pop, temp_inx, bound, 1, kk, pro, alg, rnd);
							for (size_t k = 0; k < off.size(); ++k) {
								all_off.emplace_back(off[k]);
							}
							operator_ratio[1] += 1.;
							break;
						}
					}
					if (!flag) {
						auto off = sampleInSolution(getPop()[i][j], 1, pro, alg, rnd);
						for (size_t k = 0; k < off.size(); ++k) {
							all_off.emplace_back(off[k]);
						}
						operator_ratio[0] += 1.;
					}
				}
				operator_ratio[0] /= getPop()[i].size();
				operator_ratio[1] /= getPop()[i].size();
			}
			else if (type[i] == 3) {
				//连通子空间采样
				for (size_t j = 0; j < exploit_num; ++j) {
					size_t cluster_inx = (size_t)std::floor(front_clusters.size() * rnd->uniform.next());
					//auto off = sampleInFrontNeighSpace(front_clusters[cluster_inx], bound, 1, pro, alg, rnd);
					auto off = sampleInFrontNeighSpace(front_clusters[cluster_inx], j, bound, 1, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
				operator_ratio[1] = 1.;
			}
			getOperatorRatio().emplace_back(operator_ratio);

			auto& interactive_sols = getInteractiveSols();
			for (size_t k = 0; k < interactive_sols.size(); ++k) {
				auto ind = *interactive_sols[k].back();
				getPop()[i].getOffspring()[k].variable() = ind.variable();
				getPop()[i].getOffspring()[k].objective() = ind.objective();
				getPop()[i].getOffspring()[k].setCounter(0);
			}
		}
		updateEE(0, num_exploit);
	}

	bool SPMOEA2_1::spaceSubdivision(Problem* pro, Random* rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		/*
		std::vector<size_t> front_spaces;
		for (size_t i = 0; i < num_spaces; ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_spaces.push_back(i);
			}
		}*/
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

	void SPMOEA2_1::updateFrontRegionLinkSpace(Problem* pro, Random* rnd) {
		//先找到前沿子空间
		std::vector<size_t> front_space = getFrontSpace();
		//auto link_space = clusterRegionFrontSpace(front_space);
		auto link_space = linearClusterFrontSpace(front_space, pro, rnd);//线性子空间聚类
		getFrontRegionLinkSpace() = link_space;
	}

	std::vector<std::vector<size_t>> SPMOEA2_1::linearClusterFrontSpace(std::vector<size_t>& frontspace, Problem* pro, Random* rnd) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//标记
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			//首先选择子空间内个体数最多的子空间开始聚
			std::vector<size_t> num_front_sols;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 1) {
					num_front_sols.push_back(0);
				}
				else {
					num_front_sols.push_back(getMO_HLC().getSubspaceInfo(frontspace[i]).m_front_sol_in_subspace.size());
				}
			}
			auto inx = std::distance(num_front_sols.begin(), std::max_element(num_front_sols.begin(), num_front_sols.end()));
			head_cluster.push_back(frontspace[inx]);//先加入一个子空间
			select_flag[inx] = 1;

			for (size_t i = 0; i < select_flag.size(); ++i) {
				count += select_flag[i];
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
			if (getMO_HLC().getSubspaceInfo(head_cluster[0]).m_linear_flag == 0) {
				clustered.emplace_back(head_cluster);
			}
			else {
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
									//判断是否满足邻域子空间线性关系
									if (ifContinuousSpace(inx, frontspace[k])) {
										head_cluster.push_back(frontspace[k]);
										temp.push_back(frontspace[k]);
										select_flag[k] = 1;
										count++;
									}
								}
							}
						}
					}
					if (temp.empty()) {//没有新的邻域加入
						break;
					}
					else {
						temp_cluster = temp;
					}
				}
				clustered.emplace_back(head_cluster);
			}
		}
		return clustered;
	}

	void SPMOEA2_1::clusterSubspace() {
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

	std::vector<std::vector<size_t>> SPMOEA2_1::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA2_1::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA2_1::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA2_1::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA2_1::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA2_1::findSplitDim(int inx, Problem* pro) {
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

	void SPMOEA2_1::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA2_1::recordMetrics(Problem* pro, Algorithm* alg) {
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

	void SPMOEA2_1::initiObjSpace(Problem* pro) {

	}

	std::vector<std::vector<Real>> SPMOEA2_1::qLearningNonStationaryMDP(const NonStationaryMDP& mdp, Random* rnd,
		Real learning_rate = 0.1,
		Real discount_factor = 0.9,
		Real exploration_prob = 0.2,
		int num_episodes = 1000) {
		int num_states = mdp.getTransitionProbabilities()[0].size();
		int num_actions = mdp.getTransitionProbabilities()[0][0].size();
		int horizon = mdp.getTransitionProbabilities().size();

		std::vector<std::vector<Real>> Q(num_states, std::vector<Real>(num_actions, 0.0));

		for (int episode = 0; episode < num_episodes; ++episode) {
			int state = rand() % num_states;

			for (int t = 0; t < horizon; ++t) {
				// 选择动作
				int action;
				if (rand() / static_cast<Real>(RAND_MAX) < exploration_prob) {
					action = rand() % num_actions;
				}
				else {
					auto it = std::max_element(Q[state].begin(), Q[state].end());
					action = std::distance(Q[state].begin(), it);
				}

				// 获取下一个状态和奖励
				const auto& next_state_probabilities = mdp.getTransitionProbabilities()[t][state][action];
				int next_state = generateRandomChoice(next_state_probabilities, rnd);
				Real reward = mdp.getRewardFunction()[t][state][action];

				// 更新Q值
				Q[state][action] = (1 - learning_rate) * Q[state][action] +
					learning_rate * (reward + discount_factor * *std::max_element(Q[next_state].begin(), Q[next_state].end()));

				// 更新当前状态
				state = next_state;
			}
		}

		return Q;
	}

	int SPMOEA2_1::generateRandomChoice(const std::vector<Real>& probabilities, Random* rnd) {
		Real rand_num = rnd->uniform.next();
		Real cumulative_prob = 0.0;
		for (int i = 0; i < probabilities.size(); ++i) {
			cumulative_prob += probabilities[i];
			if (rand_num <= cumulative_prob) {
				return i;
			}
		}
		return probabilities.size() - 1;
	}

	void NonStationaryMDP::generateNonStationaryTransitions() {
		// 生成随时间变化的状态转移概率
		transition_probabilities.resize(horizon, std::vector<std::vector<std::vector<Real>>>(num_states, std::vector<std::vector<Real>>(num_actions, std::vector<Real>(num_states, 0.0))));
		for (int t = 0; t < horizon; ++t) {
			for (int s = 0; s < num_states; ++s) {
				for (int a = 0; a < num_actions; ++a) {
					// 在这里可以添加对状态转移概率的具体变化规律
					transition_probabilities[t][s][a] = generateDirichletDistribution(num_states);
				}
			}
		}
	}

	void NonStationaryMDP::generateNonStationaryRewards() {
		// 生成随时间变化的即时奖励
		reward_function.resize(horizon, std::vector<std::vector<Real>>(num_states, std::vector<Real>(num_actions, 0.0)));
		for (int t = 0; t < horizon; ++t) {
			for (int s = 0; s < num_states; ++s) {
				for (int a = 0; a < num_actions; ++a) {
					// 在这里可以添加对奖励函数的具体变化规律
					reward_function[t][s][a] = generateNormalDistribution(0, 1);
				}
			}
		}
	}

	std::vector<Real> NonStationaryMDP::generateDirichletDistribution(int size) {
		std::vector<Real> distribution(size, 0.0);
		Real sum = 0.0;
		for (int i = 0; i < size - 1; ++i) {
			distribution[i] = generateUniformDistribution(0, 1 - sum);
			sum += distribution[i];
		}
		distribution[size - 1] = 1 - sum;
		return distribution;
	}

	Real NonStationaryMDP::generateUniformDistribution(Real min, Real max) {
		std::uniform_real_distribution<Real> distribution(min, max);
		static std::default_random_engine generator;
		return distribution(generator);
	}

	Real NonStationaryMDP::generateNormalDistribution(Real mean, Real stddev) {
		std::normal_distribution<Real> distribution(mean, stddev);
		static std::default_random_engine generator;
		return distribution(generator);
	}
}