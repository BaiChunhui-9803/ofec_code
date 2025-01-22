#include "spmoea1_10.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_10::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_10::run_() {
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

	void SPMOEA1_10::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_10::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(2 * getPop().size() + 1);//第二个为子代
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j]);
				}
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i + 1].push_back(&getPop()[i].getOffspring()[j]);
				}
			}
			auto& his_sols = getHisSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(his_sols[i].get());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA1_10::initiVarSpace(Problem *pro) {
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

	void SPMOEA1_10::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
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
				setArchiveNum(200);
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
		//根据子空间划分聚类
		clusterSubspace();
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

	int SPMOEA1_10::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		/********************************************************************************
							   根据子种群前沿子空间的体积分配搜索资源
		********************************************************************************/
		//根据算子类型和子连通区域长度分配计算资源
		getInteractiveSols().clear();
		std::vector<size_t> assign_pop_resource;
		size_t switch_period = 2;
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
				if (m_stage_last_time > 0) {
					interactive_type.push_back(5);
				}
				else {
					interactive_type.push_back(4);
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
				if (m_stage_last_time > 5 && m_divide_granularity < 50) {
					m_divide_granularity++;
					m_stage_last_time = 0;
				}
				//spaceSubdivision(pro, rnd);
			}
			else {
				m_stage_last_time = 0;
			}
			updateFrontSpace();
			updateFrontRegionLinkSpace(pro,rnd);
			setFrontLastGens(1);
			clusterSubspace();
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
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		SPMOEA::updateArchive(archiveNum(), pro);
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


	void SPMOEA1_10::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA1_10::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_10::updateSubPopSpace(size_t pop_inx, Problem *pro, Random *rnd) {
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

	void SPMOEA1_10::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem *pro) {
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

	void SPMOEA1_10::generateOffspring(Problem *pro, Algorithm *alg, Random *rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
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
			else if (type[i] == 2) {//前沿子空间连通集交互，连通集中采样
				//亲代所在子空间采样
				std::vector<size_t> parent_spaces;
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					auto sol = getPop()[i][j].variable().vect();
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					if (parent_spaces.empty() || std::find(parent_spaces.begin(), parent_spaces.end(), space_inx) == parent_spaces.end()) {
						parent_spaces.push_back(space_inx);
					}
				}
				std::map<size_t, std::vector<size_t>> not_front_inds;//key为子空间索引，value为个体索引
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					auto sol = getPop()[i][j].variable().vect();
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					if (std::find(front_spaces.begin(), front_spaces.end(), space_inx) == front_spaces.end()) {
						if (not_front_inds[space_inx].empty()) {
							not_front_inds.insert(std::make_pair(space_inx, j));
						}
						not_front_inds[space_inx].push_back(j);
					}
				}
				size_t count = 0;
				std::vector<size_t> ind_spaces;
				for (auto jj : not_front_inds) {
					ind_spaces.push_back(jj.first);
				}
				count = ind_spaces.size();
				//种群子空间采样
				for (size_t j = 0; j < ind_spaces.size(); ++j) {
					auto off1 = sampleInSpace(ind_spaces[j], 1, pro, alg, rnd);
					for (size_t k = 0; k < off1.size(); ++k) {
						all_off.emplace_back(off1[k]);
					}
				}
				//前沿子空间采样
				size_t assign2 = (size_t)std::floor((getPop()[i].size() - (Real)count) / 2.);
				//auto off2 = sampleInLinkSpaces(front_spaces, bound, assign2, pro, rnd);
				for (size_t j = 0; j < assign2; ++j) {
					size_t sp_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
					auto off1 = sampleInSpace(front_spaces[sp_inx], 1, pro, alg, rnd);
					for (size_t k = 0; k < off1.size(); ++k) {
						all_off.emplace_back(off1[k]);
					}
				}
				//亲代交互
				size_t assign3 = getPop()[i].size() - count - assign2;
				size_t kk = 2;
				auto off3 = sampleByDE(getPop()[i], bound, assign3, kk, pro, alg, rnd);
				//为子代赋值
				/*for (size_t k = 0; k < off2.size(); ++k) {
					all_off.emplace_back(off2[k]);
				}*/
				for (size_t k = 0; k < off3.size(); ++k) {
					all_off.emplace_back(off3[k]);
				}
			}
			else if (type[i] == 3) {//连通子空间采样
				for (size_t j = 0; j < front_clusters.size(); ++j) {
					auto off = sampleInFrontNeighSpace(front_clusters[j], bound, 1 * front_clusters[j].size(), pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else if (type[i] == 4) {
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
							auto off = sampleByDE(getPop()[i], j, bound, 1, kk, pro, alg, rnd);
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
				/*for (size_t j = 0; j < no_front_inds.size(); ++j) {
					auto off = sampleInSolution(getPop()[i][no_front_inds[j]], 1, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
					ind_order.push_back(no_front_inds[j]);
				}
				operator_ratio[0] = (Real)no_front_inds.size() / getPop()[i].size();
				size_t kk = 2;
				for (size_t j = 0; j < front_inds.size(); ++j) {
					auto off = sampleByDE(getPop()[i], front_inds[j], bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
					ind_order.push_back(front_inds[j]);
				}
				operator_ratio[1] = (Real)front_inds.size() / getPop()[i].size();*/
			}
			else if (type[i] == 5) {
				////前沿子空间采样
				//for (size_t j = 0; j < exploit_num; ++j) {
				//	size_t sp_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
				//	auto off1 = sampleInSpace(front_spaces[sp_inx], 1, pro, alg, rnd);
				//	for (size_t k = 0; k < off1.size(); ++k) {
				//		all_off.emplace_back(off1[k]);
				//	}
				//}
				//连通子空间采样
				for (size_t j = 0; j < exploit_num; ++j) {
					size_t cluster_inx = (size_t)std::floor(front_clusters.size() * rnd->uniform.next());
					auto off = sampleInFrontNeighSpace(front_clusters[cluster_inx], bound, 1, pro, alg, rnd);
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
				getPop()[i].getOffspring()[k].setTimeEvaluate(0);
			}
		}
		updateEE(0, num_exploit);
	}

	bool SPMOEA1_10::spaceSubdivision(Problem *pro, Random *rnd) {
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

	void SPMOEA1_10::clusterSubspace() {
		//每次聚类只加入一层，允许类间有重叠
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		std::vector<size_t> cluster;
		std::vector<std::vector<size_t>> clusters;
		for (size_t i = 0; i < space_num; ++i) {
			cluster.push_back(i);
		}
		clusters.emplace_back(cluster);
		getMO_HLC().getClusters().emplace_back(clusters);
	}

	std::vector<std::vector<size_t>> SPMOEA1_10::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_10::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA1_10::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_10::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem *pro, Random *rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA1_10::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_10::findSplitDim(int inx, Problem *pro) {
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

	void SPMOEA1_10::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_10::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		//使用archive计算性能指标
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			temp_pop.append(*m_archive[i]);
		}
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			pop_objs.push_back(m_archive[i]->objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);
		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

	void SPMOEA1_10::initiObjSpace(Problem *pro) {

	}
}