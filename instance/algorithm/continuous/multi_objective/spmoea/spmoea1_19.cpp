#include "spmoea1_19.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_19::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_19::run_() {
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

	void SPMOEA1_19::record() {
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
	void SPMOEA1_19::updateBuffer() {
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

	void SPMOEA1_19::initiVarSpace(Problem* pro) {
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

	void SPMOEA1_19::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
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
		m_divide_granularity = getSplitGranularity();//小于等于5
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

	int SPMOEA1_19::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		/********************************************************************************
							   根据前沿子空间的分布特点分配搜索资源演化
		********************************************************************************/
		//根据算子类型和子连通区域长度分配计算资源
		getInteractiveSols().clear();
		size_t switch_period = 5;
		//设置探索与开发的比例，由前沿子空间与非前沿子空间的个数或是体积比例来控制
		std::vector<size_t> assign_pop_resource;
		PopResourceAssign(assign_pop_resource, switch_period, pro);
		generateOffspring(pro, alg, rnd, switch_period);
		////基于个体生成子代
		//std::vector<int> interactive_type;//可以设置不同连通子空间的交互类型
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	if ((m_divide_iteration+1) % switch_period == 0) {
		//		interactive_type.push_back(1);
		//	}
		//	else {
		//		interactive_type.push_back(2);
		//	}
		//}
		//generateIndOffspring(pro, alg, rnd, switch_period, interactive_type);

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
		updateFrontSpace();//使用子种群子代更新子空间信息
		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		auto pre_num_spaces = getMO_HLC().numSubspace();
		//分析子空间是否为单段来决定是否划分
		auto front_spaces = getFrontSpace();
		size_t switch_period2 = 10;
		if (m_divide_iteration % switch_period2 == 0) {
			for (size_t i = 0; i < front_spaces.size(); ++i) {
				Real probability = ifSingleSegSpace(front_spaces[i]);
				Real rand = rnd->uniform.next();
				if (rand > probability) {//判定为多段
					spaceSubdivision(front_spaces[i], pro, rnd);
				}
			}
			updateFrontSpace();
			//更新子空间邻域
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
			}
		}
		//更新子空间连续性
		auto cur_front_spaces = getFrontSpace();
		for (size_t i = 0; i < cur_front_spaces.size(); ++i) {
			updateLinkSubspace(cur_front_spaces[i], cur_front_spaces,rnd);
		}

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

		testCoverage(pro);

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
		//两种选择模式：1、基于个体迭代的选择；2、
		//采用局部比较的方式，选择个体
		//在进行目标空间排序选择时，考虑解在决策空间的距离
		//multiObjSelection(pro);
		//sparseSelection(pro);
		//localSelection(pro);
		//localCrowdSelection(pro,rnd);
		//localCrowdSelection2(pro, rnd);
		//ensembleSelection(getPop()[0].size(), pro, rnd);
		//在每个linkSpace内结合历史非支配解选择保留的个体
		if ((m_divide_iteration + 1) % switch_period != 0) {
			SolutionSelection(pro, rnd);
		}
		

		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		//SPMOEA::recordMetrics(pro, alg);
		SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA1_19::testCoverage(Problem* pro) {
		//检测子空间是否包含不同的PS片段
		m_real_front_subspaces.clear();
		m_multi_segment_subspaces.clear();
		m_match_subspace.clear();
		m_error_subspace.clear();
		m_loss_subspace.clear();
		m_match_multi_seg_subspace.clear();
		auto front_sp = getFrontSpace();
		//真实PS所在的子空间
		//std::vector<size_t> real_front_subspaces;
		////具有多段的子空间数目
		//std::vector<size_t> multi_segment_subspaces;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			auto& box = getMO_HLC().subspaceTree().getBox(i);
			//沿着第一维构造解
			size_t num = 100;
			std::vector<size_t> flag(num - 1, 1);
			Real delta = (box[0].second - box[0].first) / num;
			for (size_t j = 0; j < num - 1; ++j) {
				std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
				auto sol = CAST_CONOP(pro)->createVar(s);

				//判断sol是否在box内
				for (size_t k = 0; k < box.size(); ++k) {
					/*if (dynamic_cast<const VariableVector<>&>(sol->variableBase())[k] > box[k].second || dynamic_cast<const VariableVector<>&>(sol->variableBase())[k] < box[k].first) {
						flag[j] = 0;
						break;
					}*/
					if (sol[k] > box[k].second || sol[k] < box[k].first) {
						flag[j] = 0;
						break;
					}
				}
				/*std::vector<Real> sol(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
				for (size_t k = 0; k < box.size(); ++k) {
					if (sol[k] > box[k].second || sol[k] < box[k].first) {
						flag[j] = 0;
						break;
					}
				}*/
			}
			//统计段数
			size_t segment_time = 0;
			for (size_t j = 1; j < flag.size(); ++j) {
				if (flag[0] == 0) {
					if (flag[j - 1] == 0 && flag[j] == 1) {
						segment_time++;
					}
				}
				else {
					if (flag[j - 1] == 1 && flag[j] == 0) {
						segment_time++;
					}
				}
			}
			if (segment_time > 0) {
				m_real_front_subspaces.push_back(i);
				if (segment_time > 1) {
					m_multi_segment_subspaces.push_back(i);
				}
			}
			//std::cout << "the number of segment of the " << i << " th subspace is: " << segment_time << std::endl;
			//1、看子空间的非线性程度

			//2、看子空间内的解是否可分为多个类，采用dbscan聚类

		}
		//真实子空间与逼近子空间的差异
		/*std::vector<size_t> match_subspace;
		std::vector<size_t> error_subspace;
		std::vector<size_t> loss_subspace;
		std::vector<size_t> match_multi_seg_subspace;*/
		for (size_t i = 0; i < front_sp.size(); ++i) {
			if (std::find(m_real_front_subspaces.begin(), m_real_front_subspaces.end(), front_sp[i]) != m_real_front_subspaces.end()) {
				m_match_subspace.push_back(front_sp[i]);
			}
			else {
				m_error_subspace.push_back(front_sp[i]);
			}
		}
		for (size_t i = 0; i < m_real_front_subspaces.size(); ++i) {
			if (std::find(front_sp.begin(), front_sp.end(), m_real_front_subspaces[i]) == front_sp.end()) {
				m_loss_subspace.push_back(m_real_front_subspaces[i]);
			}
		}
		for (size_t i = 0; i < m_match_subspace.size(); ++i) {
			if (std::find(m_multi_segment_subspaces.begin(), m_multi_segment_subspaces.end(), m_match_subspace[i]) != m_multi_segment_subspaces.end()) {
				m_match_multi_seg_subspace.push_back(m_match_subspace[i]);
			}
		}
		std::cout << "the number of total   subspace is: " << getMO_HLC().numSubspace() << std::endl;
		std::cout << "the number of real-PS subspace is: " << m_real_front_subspaces.size() << std::endl;
		std::cout << "the number of mul-seg subspace is: " << m_multi_segment_subspaces.size() << std::endl;
		std::cout << "the number of matched subspace is: " << m_match_subspace.size() << std::endl;
		std::cout << "match multiseg subspace number is: " << m_match_multi_seg_subspace.size() << std::endl;
		std::cout << "the number of errored subspace is: " << m_error_subspace.size() << std::endl;
		std::cout << "the number of lossed  subspace is: " << m_loss_subspace.size() << std::endl;

	}

	void SPMOEA1_19::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA1_19::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_19::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
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

	void SPMOEA1_19::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		std::vector<Real> ratio_ee;
		auto front_spaces = getFrontSpace();
		std::vector<size_t> behind_space;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			if (std::find(front_spaces.begin(), front_spaces.end(), j) == front_spaces.end()) {
				behind_space.push_back(j);
			}
		}
		ratio_ee.push_back((Real)behind_space.size()/ getMO_HLC().numSubspace());
		ratio_ee.push_back(1.-(Real)behind_space.size() / getMO_HLC().numSubspace());
		m_num_e_e.emplace_back(ratio_ee);
		//updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_19::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, size_t switch_period) {
		size_t num_exploit = 0;
		auto bound = CAST_CONOP(pro)->boundary();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		for (size_t i = 0; i < getPop().size(); ++i) {
			//计算资源分配，根据种群的历史前沿解在边界内的比例分配计算资源,还有非前沿子空间
			size_t explore_num = getPopsize()*getEERatio().back()[0];
			size_t exploit_num = getPopsize() - explore_num;
			auto front_spaces = getFrontSpace();
			std::vector<size_t> behind_space;
			for (size_t j = 0; j < getPop().size(); ++j) {
				for (size_t k = 0; k < getPop()[j].size(); ++k) {
					auto& sol = getPop()[j][k].variable().vect();
					auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					if (std::find(front_spaces.begin(), front_spaces.end(), inx) == front_spaces.end()) {
						behind_space.push_back(inx);
					}
				}
			}
			std::vector<size_t> temp_num;
			if ((m_divide_iteration + 1) % switch_period == 0) {
				//模型采样，根据子空间内前沿解的稀疏程度分配采样个数


			}
			else {
				//演化进化，根据子空间的连续集分配计算资源
				//考察每个子空间连续集中当前的个体数，按比例分配，最小为1
				std::vector<size_t> space_pop_num;
				size_t total_num = 0;
				auto& his_front_sol = getHisFrontSols();
				for (size_t i = 0; i < front_spaces.size(); ++i) {
					std::vector<size_t> linkspaces = getMO_HLC().getSubspaceInfo(front_spaces[i]).m_link_subspaces;
					linkspaces.push_back(front_spaces[i]);
					size_t count = 0;
					//历史非支配解在连通集中的个数
					for (size_t j = 0; j < his_front_sol.size(); ++j) {
						auto& sol = his_front_sol[j]->variable().vect();
						auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
						if (std::find(linkspaces.begin(), linkspaces.end(), inx) != linkspaces.end()) {
							count++;
						}
					}
					total_num += count;
					space_pop_num.push_back(count);
				}
				//for (size_t i = 0; i < front_spaces.size(); ++i) {
				//	std::vector<size_t> linkspaces = getMO_HLC().getSubspaceInfo(front_spaces[i]).m_link_subspaces;
				//	linkspaces.push_back(front_spaces[i]);
				//	size_t count = 0;
				//	//历史非支配解在连通集中的个数
				//	for (size_t j = 0; j < getPop().size(); ++j) {
				//		for (size_t k = 0; k < getPop()[j].size(); ++k) {
				//			auto& sol = getPop()[j][k].variable().vect();
				//			auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				//			if (std::find(linkspaces.begin(), linkspaces.end(), inx) != linkspaces.end()) {
				//				count++;
				//			}
				//		}
				//	}
				//	total_num += count;
				//	space_pop_num.push_back(count);
				//}
				size_t total = 0;
				for (size_t i = 0; i < front_spaces.size(); ++i) {
					size_t assign_num = (size_t)std::ceil(exploit_num * space_pop_num[i] / (Real)total_num);
					total += assign_num;
					temp_num.push_back(assign_num);
				}
				while (total > exploit_num) {
					size_t inx = std::distance(temp_num.begin(), std::max_element(temp_num.begin(), temp_num.end()));
					if (temp_num[inx] == 1) {
						break;
					}
					else {
						temp_num[inx]--;
						total--;
					}
				}
			}
			size_t total_pop_num = 0;
			for (size_t j = 0; j < temp_num.size(); ++j) {
				total_pop_num += temp_num[j];
			}
			explore_num = getPopsize() - total_pop_num;
			//子代个体数量处理
			while (getPop()[i].getOffspring().size() - getPop()[i].size() > getPopsize()) {
				getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
			}
			while (getPop()[i].getOffspring().size() - getPop()[i].size() < getPopsize()) {
				Solution<> temp_ind(M, CAST_CONOP(pro)->numberConstraints(), bound.size());
				getPop()[i].getOffspring().append(temp_ind);
			}
			//auto bound = CAST_CONOP(pro)->boundary();
			auto front_clusters = getFrontRegionLinkSpace();
			
			size_t total_spaces = getMO_HLC().numSubspace();
			std::vector<size_t> space_fre;
			for (size_t j = 0; j < total_spaces; ++j) {
				space_fre.push_back(getMO_HLC().getSubspaceInfo(j).m_sub_freq);
			}
			std::vector<Real> operator_ratio(2, 0.);
			std::vector<std::vector<Real>> all_off;

			//根据分配的探索与开发比例生成子代
			if ((m_divide_iteration + 1) % switch_period == 0) {
				//模型采样，根据子空间内前沿解的稀疏程度分配采样个数


			}
			else {//在子空间linkspace间交互,linkspace的前沿解构成子种群交互
				//对于前沿子空间，exploit
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					std::vector<size_t> linkspaces = getMO_HLC().getSubspaceInfo(front_spaces[j]).m_link_subspaces;
					linkspaces.push_back(front_spaces[j]);
					auto off = sampleInFrontNeighSpace(linkspaces, bound, temp_num[j], pro, alg, rnd);
					//auto off = sampleInFrontNeighSpace(linkspaces, j, bound, pop_resource[j], pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
				//对于非前沿子空间，子空间内部交互，解不足，变异，解充足，DE
				//依次选择子空间内个体数最少的子空间生成新的解
				std::vector<size_t> sol_num;
				for (size_t j = 0; j < behind_space.size(); ++j) {
					sol_num.push_back(getMO_HLC().getSubspaceInfo(behind_space[j]).m_history_inds.size());
				}
				for (size_t j = 0; j < explore_num; ++j) {
					size_t sele_inx = std::distance(sol_num.begin(), std::min_element(sol_num.begin(), sol_num.end()));
					//在子空间中根据历史解个数生成解
					size_t space_sol_num = getMO_HLC().getSubspaceInfo(behind_space[sele_inx]).m_history_inds.size();
					if (space_sol_num < 5) {
						//子空间内随机采样
						auto off = sampleInSubspace(behind_space[sele_inx],1,pro,alg,rnd);
						for (size_t k = 0; k < off.size(); ++k) {
							all_off.emplace_back(off[k]);
						}
					}
					else {
						//子空间前沿解
						SPMOEA_pop temp_pop(0, pro);
						auto& his_sols = getMO_HLC().getSubspaceInfo(sele_inx).m_history_inds;
						for (size_t k = 0; k < his_sols.size();++k) {
							temp_pop.append(*his_sols[k]);
						}
						//在合成种群中演化
						size_t kk = 2;
						auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
						for (size_t k = 0; k < off.size(); ++k) {
							all_off.emplace_back(off[k]);
						}
					}
					sol_num[sele_inx]++;
				}
			}
		    //operator_ratio[0] /= getPop()[i].size();
			//operator_ratio[1] /= getPop()[i].size();
			//getOperatorRatio().emplace_back(operator_ratio);

			auto& interactive_sols = getInteractiveSols();
			for (size_t k = 0; k < interactive_sols.size(); ++k) {
				auto ind = *interactive_sols[k].back();
				getPop()[i].getOffspring()[k].variable() = ind.variable();
				getPop()[i].getOffspring()[k].objective() = ind.objective();
				getPop()[i].getOffspring()[k].setTimeEvaluate(0);
			}
		}
		//updateEE(0, num_exploit);
	}

	void SPMOEA1_19::generateIndOffspring(Problem* pro, Algorithm* alg, Random* rnd, size_t switch_peroid, std::vector<int> type) {
		size_t num_exploit = 0;
		auto bound = CAST_CONOP(pro)->boundary();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		//前面的种群在自己的连通集之间扩展
		for (size_t i = 0; i < getPop().size(); ++i) {
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
			if (type[i] == 1) {//模型采样
				//模型采样，根据子空间内前沿解的稀疏程度分配采样个数

			}
			else if (type[i] == 2) {//根据是否为前沿个体进行采样
				//基于每个个体采样
				size_t kk = 2;
				auto& his_front_sols = getHisFrontSols();
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					bool flag = false;
					auto& sol = getPop()[i][j].variable().vect();
					size_t space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					for (size_t k = 0; k < his_front_sols.size(); ++k) {
						auto& sol2 = his_front_sols[k]->variable().vect();
						if (ifSame(sol, sol2)) {
							flag = true;
							break;
						}
					}
					if (flag) {//基于连通集采样
						std::vector<size_t> linkspaces = getMO_HLC().getSubspaceInfo(space_inx).m_link_subspaces;
						linkspaces.push_back(space_inx);
						auto off = sampleInFrontNeighSpace(linkspaces, bound, 1, pro, alg, rnd);
						//auto off = sampleInFrontNeighSpace(linkspaces, j, bound, pop_resource[j], pro, alg, rnd);
						for (size_t k = 0; k < off.size(); ++k) {
							all_off.emplace_back(off[k]);
						}
						operator_ratio[0]= operator_ratio[0]+1.;
					}
					else {//与邻域子空间前沿解交互
						auto& neighs = getMO_HLC().getSubspaceInfo(space_inx).m_sub_neighbors;
						std::vector<size_t> neis;
						for (auto jj : neighs) {
							neis.push_back(jj);
						}
						//neis.push_back(space_inx);
						//组成新种群，基于个体DE差分
						SPMOEA_pop temp_pop(0, pro);
						for (size_t k = 0; k < neis.size(); ++k) {
							auto& sp_front_sols = getMO_HLC().getSubspaceInfo(neis[k]).m_subspace_front_sol;
							for (size_t p = 0; p < sp_front_sols.size(); ++p) {
								temp_pop.append(*sp_front_sols[p]);
							}
						}
						auto& sp_front_sols = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
						for (size_t p = 0; p < sp_front_sols.size(); ++p) {
							auto& sol2 = sp_front_sols[p]->variable().vect();
							if (!ifSame(sol, sol2)) {
								temp_pop.append(*sp_front_sols[p]);
							}
						}
						temp_pop.append(getPop()[i][j]);
						auto off = sampleByDE(temp_pop, temp_pop.size() - 1, bound, 1, kk, pro, alg, rnd);
						for (size_t k = 0; k < off.size(); ++k) {
							all_off.emplace_back(off[k]);
						}
						operator_ratio[1]= operator_ratio[1] + 1.;
					}
				}
				operator_ratio[0] /= (Real)getPop()[i].size();
				operator_ratio[1] /= (Real)getPop()[i].size();
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

	bool SPMOEA1_19::spaceSubdivision(size_t space_inx, Problem* pro, Random* rnd) {
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
			auto divide_info = findSplitDimAndPos(space_inx, pro);
			splitSpace(space_inx, 2, divide_info.first, divide_info.second, false, pro, rnd);
		}
		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}
		return flag;
	}

	std::vector<std::vector<size_t>> SPMOEA1_19::linearClusterFrontSpace(std::vector<size_t>& frontspace, Problem* pro, Random* rnd) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//标记
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space=frontspace[i];
					select_flag[i] = 1;
				}
			}
			head_cluster.push_back(begin_space);//先加入一个子空间
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
					auto& linkspace = getMO_HLC().getSubspaceInfo(inx).m_link_subspaces;
					for (size_t k = 0; k < frontspace.size(); ++k) {
						if (select_flag[k] == 0) {
							if (std::find(linkspace.begin(), linkspace.end(), frontspace[k]) != linkspace.end()) {
								head_cluster.push_back(frontspace[k]);
								temp.push_back(frontspace[k]);
								select_flag[k] = 1;
								count++;
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
		return clustered;
	}

	void SPMOEA1_19::updateFrontRegionLinkSpace(Problem* pro, Random* rnd) {
		//先找到前沿子空间
		std::vector<size_t> front_space = getFrontSpace();
		//auto link_space = clusterRegionFrontSpace(front_space);
		auto link_space = linearClusterFrontSpace(front_space,pro,rnd);//线性子空间聚类
		getFrontRegionLinkSpace() = link_space;
	}

	void SPMOEA1_19::updateLinkSubspace(size_t inx,const std::vector<size_t>& front_spaces,Random* rnd) {
		auto& neigh = getMO_HLC().getSubspaceInfo(inx).m_sub_neighbors;
		std::vector<size_t> link_neighs;
		for (auto& i : neigh) {
			if (std::find(front_spaces.begin(), front_spaces.end(), i) != front_spaces.end()) {
				Real link_probability = ifContinuousSpace(inx,i);
				Real rand = rnd->uniform.next();
				if (rand <= link_probability) {
					link_neighs.push_back(i);
				}
			}
		}
		getMO_HLC().getSubspaceInfo(inx).m_link_subspaces = link_neighs;
	}

	void SPMOEA1_19::clusterSubspace() {
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

	std::vector<std::vector<size_t>> SPMOEA1_19::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_19::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA1_19::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_19::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA1_19::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_19::findSplitDim(int inx, Problem* pro) {
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

	void SPMOEA1_19::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_19::recordMetrics(Problem* pro, Algorithm* alg) {
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

	void SPMOEA1_19::initiObjSpace(Problem* pro) {

	}
}