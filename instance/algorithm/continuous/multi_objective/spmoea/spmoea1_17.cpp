#include "spmoea1_17.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_17::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		auto& v = *m_param;
		
		m_num_push = v.get<int>("number of push");
		m_num_extend = v.get<int>("number of extend");
		m_switch_period = m_num_push+m_num_extend;
		//m_neighs = 5;
	}

	void SPMOEA1_17::run_() {
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

	void SPMOEA1_17::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_17::updateBuffer() {
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
			auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(his_sols[i].get());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA1_17::initiVarSpace(Problem *pro) {
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

	void SPMOEA1_17::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
		//在各个子空间生成一个子种群
		size_t number_objectives = CAST_CONOP(pro)->numberObjectives();
		if (number_objectives == 2) {
			setArchiveNum(100);
		}
		else if (number_objectives == 3) {
			setArchiveNum(300);
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
		auto num_var = CAST_CONOP(pro)->numberVariables();
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
		updateVarSpaceRank(pro,rnd);
		updateEE(pop_num, 0);
		SPMOEA::recordMetrics(pro, alg);

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
		//m_divide_granularity = getSplitGranularity();
		m_divide_granularity = 20 - 10. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
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

		//真实PS与子空间之间的关系
		//沿着第一维直接构造解
		size_t num = 10000;
		//auto bound = CAST_CONOP(pro)->boundary();
		Real delta = (bound[0].second - 0.001 - bound[0].first) / num;
		std::vector<size_t> real_front;
		std::vector<size_t> multiseg_front;
		size_t temp_space = INT16_MAX;
		for (size_t j = 0; j < num - 1; ++j) {
			std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), bound[0].first + 0.0001 + (j + 1) * delta);
			auto sol = CAST_CONOP(pro)->createVar(s);
			//判断sol在哪个子空间
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			if (space_inx != temp_space) {
				if (std::find(real_front.begin(), real_front.end(), space_inx) == real_front.end()) {
					real_front.push_back(space_inx);
				}
				else {
					if (space_inx != real_front.back()) {
						if (std::find(multiseg_front.begin(), multiseg_front.end(), space_inx) == multiseg_front.end()) {
							multiseg_front.push_back(space_inx);
						}
					}
				}
				m_real_all_front_subspaces.push_back(space_inx);
			}
			temp_space = space_inx;
		}
		m_real_front_subspaces = real_front;
		m_multi_segment_subspaces = multiseg_front;

		//SPMOEA::record();
	}

	int SPMOEA1_17::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		/********************************************************************************
							   根据子种群前沿子空间的体积分配搜索资源
		********************************************************************************/
		//根据算子类型和子连通区域长度分配计算资源
		getInteractiveSols().clear();
		std::vector<size_t> assign_pop_resource;
		//size_t switch_period = 3;
		PopResourceAssign(assign_pop_resource, m_switch_period, pro);
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
			//auto front_clusters = getFrontRegionLinkSpace();
			//auto front_spaces = getFrontSpace();
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (m_divide_iteration % m_switch_period < m_num_push) {
					interactive_type.push_back(4);
					//根据子空间搜索程度和子控空间排序布局新的爬山位置
					if (m_divide_iteration>0) {
						assignPop(pro, alg, rnd);
					}
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
			for (size_t j = 0; j <getPop()[i].size(); j++) {
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}
		SPMOEA::NDSort(offspring_pop);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		updateSubspaceFrontSol(offspring_pop, pro, rnd);//子代更新子空间前沿解和代表解
		updateFrontSpace();//使用子种群子代更新子空间信息
		auto pre_front_space = getFrontSpace();
		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		auto pre_num_spaces = getMO_HLC().numSubspace();
		Real total_volume = getVarSpaceVolume();
		size_t partition_period = (size_t)std::floor(4*std::sqrt(pro->numberVariables()));
		partition_period = 50000;
		//子空间细分方式：1、粒度划分；2、多段划分；3、非线性程度划分
		// 1、粒度划分
		if ((m_divide_iteration+1) % partition_period == 0) {
			auto divide_flag = spaceSubdivision(pro, rnd);
			if (divide_flag) {
				updateFrontSpace();
				//更新子空间邻域
				for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
					getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
				}
			}
		}

		//2.1、 切分多段子空间
		//for (size_t i = 0; i < pre_front_space.size(); ++i) {
		//	Real single_seg_pro = ifSingleSegSpace(pre_front_space[i]);
		//	auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(pre_front_space[i]);
		//	size_t min_num = 10;
		//	if (single_seg_pro < 0.3) {//多段，切分
		//		if (space_volume / total_volume > std::pow(1. / (Real)min_num, CAST_CONOP(pro)->numberVariables())) {
		//			int dim = findSplitDim(pre_front_space[i], pro);
		//			auto& space_bound = getMO_HLC().subspaceTree().getBox(pre_front_space[i]);
		//			Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
		//			splitSpace(pre_front_space[i], 2, dim, pos, false, pro, rnd);
		//		}
		//		
		//	}
		//}
		//2.2、切分多段子空间
		/*for (size_t i = 0; i < pre_front_space.size(); ++i) {
			auto space_status = ifMultiSegSpace(pre_front_space[i]);
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(pre_front_space[i]);
			if (space_status.first) {
				size_t dim_inx = space_status.second.first;
				Real offset = space_status.second.second;
				auto& space_bound = getMO_HLC().subspaceTree().getBox(pre_front_space[i]);
				Real pos = space_bound[dim_inx].first + offset;
				splitSpace(pre_front_space[i], 2, dim_inx, pos, false, pro, rnd);
			}
		}*/
		//更新子空间连续性
		auto cur_num_spaces = getMO_HLC().numSubspace();
		if (cur_num_spaces > pre_num_spaces) {
			updateFrontSpace();
			//更新子空间邻域
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
			}
		}
		auto cur_front_spaces = getFrontSpace();
		/*for (size_t i = 0; i < cur_front_spaces.size(); ++i) {
			updateLinkSubspace(cur_front_spaces[i], cur_front_spaces, rnd);
		}*/
		setAddNumSpace(cur_num_spaces - pre_num_spaces);
		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		if ((m_divide_iteration + 1) % 1 == 0) {
			updateVarSpaceRank(pro, rnd);
		}
		//updateVarSpaceDominanceRank(pro, rnd);
		updateFrontRegionLinkSpace(pro, rnd);
		setFrontLastGens(1);
		clusterSubspace();
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		//SPMOEA::updateArchive(archiveNum(), pro);
		//updateObjSpace();
		testCoverage(pro);

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
		//ensembleSelection(getPop()[0].size(), pro, rnd);
		if (m_divide_iteration % m_switch_period < m_num_push) {
			SolutionSelection(pro, rnd);
		}
		
		////k近邻选择
		//auto sele_inx=neighLocalSelection(5,pro);
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//		getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
		//	}
		//}

		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		//SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA1_17::assignPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//先得到非前沿子空间
		std::vector<size_t> behind_spaces;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			if (getMO_HLC().getSubspaceInfo(j).m_best_rank > 0) {
				behind_spaces.push_back(j);
			}
		}
		std::vector<Real> box_span;
		for (size_t i = 0; i < behind_spaces.size(); ++i) {
			auto box = getMO_HLC().subspaceTree().getBox(behind_spaces[i]);
			Real mean_span = 0.;
			for (size_t j = 0; j < box.size(); ++j) {
				mean_span += (box[j].second - box[j].first);
			}
			mean_span /= box.size();
			box_span.push_back(mean_span);
		}
		std::vector<Real> space_density;
		for (size_t j = 0; j < behind_spaces.size(); ++j) {
			space_density.push_back(getMO_HLC().getSubspaceInfo(behind_spaces[j]).m_sub_freq / box_span[j]);
		}
		std::vector<Real> space_rank;
		for (size_t j = 0; j < behind_spaces.size(); ++j) {
			space_rank.push_back(1 + getMO_HLC().getSubspaceInfo(behind_spaces[j]).m_best_rank);
		}
		//子空间排序
		std::vector<std::vector<Real>*> objs;
		std::vector<std::vector<Real>> all_metrics;
		for (size_t i = 0; i < space_rank.size(); ++i) {
			std::vector<Real> temp = {space_density[i],space_rank[i]};
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//第一层子控空间
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(behind_spaces[i]);
			}
		}
		//将种群中到前沿的个体初始化
		std::vector<size_t> front_space=getFrontSpace();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& sol = getPop()[i][j].variable().vect();
				auto space_index = getMO_HLC().subspaceTree().getRegionIdx(sol);
				auto rank_sp = getMO_HLC().getSubspaceInfo(space_index).m_best_rank;
				if (rank_sp == 0) {
					auto& front_sols = getMO_HLC().getSubspaceInfo(space_index).m_front_sol_in_subspace;
					for (size_t k = 0; k < front_sols.size(); ++k) {
						auto& sol2 = front_sols[k]->variable().vect();
						if (ifSame(sol, sol2)) {
							//初始化到一个子空间
							std::vector<size_t> space;
							for (size_t p = 0; p < first_layer.size(); ++p) {
								auto& space_front_sol = getMO_HLC().getSubspaceInfo(first_layer[p]).m_subspace_front_sol;
								if (space_front_sol.size() > 0) {
									space.push_back(first_layer[p]);
								}
							}
							size_t space_inx;
							if (space.size() > 0) {
								space_inx = space[(size_t)std::floor(space.size() * rnd->uniform.next())];
							}
							else if (behind_spaces.size()==0) {
								do {
									space_inx = front_space[(size_t)std::floor(front_space.size() * rnd->uniform.next())];
								} while (getMO_HLC().getSubspaceInfo(space_inx).m_history_inds.size() == 0);
							}
							else{
								do {
									space_inx = behind_spaces[(size_t)std::floor(behind_spaces.size() * rnd->uniform.next())];
								} while (getMO_HLC().getSubspaceInfo(space_inx).m_history_inds.size()==0);
							}
							auto& space_front_sol = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
							size_t ind_inx = (size_t)std::floor(space_front_sol.size() * rnd->uniform.next());
							getPop()[i][j].variable() = space_front_sol[ind_inx]->variable();
							getPop()[i][j].objective() = space_front_sol[ind_inx]->objective();

							//if (space_front_sol.size() < 1) {
							//	//std::vector<Real> new_sol;
							//	//auto& box = getMO_HLC().subspaceTree().getBox(space);
							//	//for (size_t p = 0; p < pro->numberVariables(); ++p) {
							//	//	new_sol.push_back(box[p].first+rnd->uniform.next()*(box[p].second-box[p].first));
							//	//}
							//	//Solution<> ind(pro->numberObjectives(), 0);
							//	//ind.variable().vect() = new_sol;
							//	//ind.evaluate(pro,alg);
							//	//getPop()[i][j].variable() = ind.variable();
							//	//getPop()[i][j].objective() = ind.objective();
							//	//SPMOEA_pop temp_pop(0, pro);
							//	//temp_pop.append(ind);
							//	//SPMOEA::updateHistoryInfo(temp_pop, pro);
							//	//updateSubspaceFrontSol(temp_pop, pro, rnd);//子代更新子空间前沿解和代表解
							//}
							//else {
							//	//随机选择一个子空间前沿解
							//	size_t ind_inx = (size_t)std::floor(space_front_sol.size() * rnd->uniform.next());
							//	getPop()[i][j].variable() = space_front_sol[ind_inx]->variable();
							//	getPop()[i][j].objective() = space_front_sol[ind_inx]->objective();
							//}
							m_immigrant_times++;
							break;
						}
					}
				}
			}
		}
	}

	void SPMOEA1_17::testCoverage(Problem* pro) {
		//检测子空间是否包含不同的PS片段
		//m_real_front_subspaces.clear();
		//m_multi_segment_subspaces.clear();
		m_match_subspace.clear();
		m_error_subspace.clear();
		m_loss_subspace.clear();
		m_match_multi_seg_subspace.clear();
		auto front_sp = getFrontSpace();
		////动态分割下使用下面的代码，真实PS所在的子空间
		////std::vector<size_t> real_front_subspaces;
		//////具有多段的子空间数目
		////std::vector<size_t> multi_segment_subspaces;
		//for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//	auto& box = getMO_HLC().subspaceTree().getBox(i);
		//	//沿着第一维构造解
		//	size_t num = 200;
		//	std::vector<size_t> flag(num - 1, 1);
		//	Real delta = (box[0].second - box[0].first) / num;
		//	for (size_t j = 0; j < num - 1; ++j) {
		//		std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
		//		auto sol = CAST_CONOP(pro)->createVar(s);
		//		//判断sol是否在box内
		//		for (size_t k = 0; k < box.size(); ++k) {
		//			if (sol[k] > box[k].second || sol[k] < box[k].first) {
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
		//		m_real_front_subspaces.push_back(i);
		//		if (segment_time > 1) {
		//			m_multi_segment_subspaces.push_back(i);
		//		}
		//	}
		//	//std::cout << "the number of segment of the " << i << " th subspace is: " << segment_time << std::endl;
		//	//1、看子空间的非线性程度

		//	//2、看子空间内的解是否可分为多个类，采用dbscan聚类

		//}

		//测试计算邻域识别的准确性
		if (m_divide_iteration % m_switch_period >= m_num_push) {
			std::cout << "************************************************************** " << std::endl;
			std::cout << "neighs match result: " << std::endl;
			int total_match_num = 0;
			for (size_t i = 0; i < m_recog_neighs.size(); ++i) {
				int match_num = 0;
				auto center_space = m_recog_neighs[i][0];
				if (std::find(m_real_front_subspaces.begin(), m_real_front_subspaces.end(), center_space) != m_real_front_subspaces.end()) {
					//找哪些位置匹配
					std::vector<size_t> match_inx;
					for (size_t j = 0; j < m_real_all_front_subspaces.size(); ++j) {
						if (m_real_all_front_subspaces[j] == center_space) {
							match_inx.push_back(j);
						}
					}
					for (size_t j = 0; j < match_inx.size(); ++j) {
						if (match_inx[j] == 0) {
							size_t pos1 = m_real_all_front_subspaces[match_inx[j] + 1];
							for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
								if (m_recog_neighs[i][k] == pos1) {
									match_num++;
								}
							}
						}
						else if (match_inx[j] == m_real_all_front_subspaces.size() - 1) {
							size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
							for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
								if (m_recog_neighs[i][k] == pos1) {
									match_num++;
								}
							}
						}
						else {
							size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
							size_t pos2 = m_real_all_front_subspaces[match_inx[j] + 1];
							for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
								if (m_recog_neighs[i][k] == pos1 || m_recog_neighs[i][k] == pos2) {
									match_num++;
								}
							}
						}
					}
					total_match_num += match_num;
				}
				else {
					match_num = -1;
				}
				std::cout << match_num <<" ";
			}
			std::cout <<std::endl;
			//输出总的匹配率
			std::cout << "neighs match result: " <<std::setprecision(2)<< (Real)total_match_num/(2.*m_real_front_subspaces.size()-2) << std::endl;
		}
		
		//真实子空间与逼近子空间的差异
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
		
		std::cout << "************************************************************** " << std::endl;
		std::cout << "the number of total   subspace is: " << getMO_HLC().numSubspace() << std::endl;
		std::cout << "the number of real-PS subspace is: " << m_real_front_subspaces.size() << std::endl;
		std::cout << "the number of mul-seg subspace is: " << m_multi_segment_subspaces.size() << std::endl;
		std::cout << "the number of matched subspace is: " << m_match_subspace.size() << std::endl;
		std::cout << "match multiseg subspace number is: " << m_match_multi_seg_subspace.size() << std::endl;
		std::cout << "the number of errored subspace is: " << m_error_subspace.size() << std::endl;
		std::cout << "the number of lossed  subspace is: " << m_loss_subspace.size() << std::endl;
		
		std::cout << "************************************************************** " << std::endl;
		std::cout << "space sequence: " << std::endl;
		for (size_t i = 0; i < m_real_all_front_subspaces.size(); ++i) {
			std::cout << m_real_all_front_subspaces[i] << " ";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "error index: " << std::endl;
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			std::cout << m_error_subspace[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "frequency  : ";
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			size_t fre = getMO_HLC().getSubspaceInfo(m_error_subspace[i]).m_sub_freq;
			std::cout << fre << " ";
		}
		std::cout << std::endl;
		std::cout << "rank       : ";
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			size_t rank = getMO_HLC().getSubspaceInfo(m_error_subspace[i]).m_best_rank;
			std::cout << rank << " ";
		}
		std::cout << std::endl;
		
		std::cout << "************************************************************** " << std::endl;
		std::cout << "loss index: " << std::endl;
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			std::cout << m_loss_subspace[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "frequency : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			size_t fre = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_sub_freq;
			std::cout << fre << " ";
		}
		std::cout << std::endl;
		std::cout << "rank      : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			size_t rank = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_best_rank;
			std::cout << rank << " ";
		}
		std::cout << std::endl;
		std::cout << "pop span  : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			auto &his_sols = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_history_inds;
			auto& box = getMO_HLC().subspaceTree().getBox(m_loss_subspace[i]);
			std::vector<Real> span;
			if (his_sols.size() > 0) {
				for (size_t k = 0; k < pro->numberVariables(); ++k) {
					Real min_v = (Real)INT16_MAX;
					Real max_v = -1. * INT16_MAX;
					for (size_t j = 0; j < his_sols.size(); ++j) {
						if (his_sols[j]->variable()[k] > max_v) {
							max_v = his_sols[j]->variable()[k];
						}
						if (his_sols[j]->variable()[k] < min_v) {
							min_v = his_sols[j]->variable()[k];
						}
					}
					Real temp = max_v - min_v;
					span.push_back(temp / (box[k].second - box[k].first));
				}
			}
			else {
				for (size_t k = 0; k < pro->numberVariables(); ++k) {
					span.push_back(0.);
				}
			}
			
			auto max_inx = std::distance(span.begin(), max_element(span.begin(), span.end()));
			auto min_inx = std::distance(span.begin(), min_element(span.begin(), span.end()));
			std::cout <<std::setprecision(2)<<"("<<span[min_inx] << "," << span[max_inx] << ")";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "the   number   of   immigrant  is: " << m_immigrant_times << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
	}


	void SPMOEA1_17::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA1_17::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_17::updateSubPopSpace(size_t pop_inx, Problem *pro, Random *rnd) {
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

	void SPMOEA1_17::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem *pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		//updateFrontRegionLinkSpace();
		/*if (m_stage_last_time > 0) {
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
		}*/
		assign_pop_resource.push_back(max_pop_size);
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_17::generateOffspring(Problem *pro, Algorithm *alg, Random *rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
		size_t num_exploit = 0;
		auto search_bound = CAST_CONOP(pro)->boundary();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		//前面的种群在自己的连通集之间扩展
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t exploit_num = pop_resource[i];
			num_exploit += exploit_num;
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
				m_recog_neighs.clear();
				std::vector<Real> box_span;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					auto box = getMO_HLC().subspaceTree().getBox(front_spaces[j]);
					Real mean_span = 0.;
					for (size_t j = 0; j < box.size(); ++j) {
						mean_span += (box[j].second - box[j].first);
					}
					mean_span /= box.size();
					box_span.push_back(mean_span);
				}
				std::vector<size_t> front_fre;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					front_fre.push_back(getMO_HLC().getSubspaceInfo(front_spaces[j]).m_sub_freq);
				}
				std::vector<Real> space_density;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					space_density.push_back(getMO_HLC().getSubspaceInfo(front_spaces[j]).m_sub_freq / box_span[j]);
				}
				std::vector<std::map<size_t, Real>> all_front_nei_prob;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					auto space_prob = frontSpaceLinkProbability(front_spaces[j]);
					all_front_nei_prob.emplace_back(space_prob);
					/***************测试识别的邻域关系是否准确**************/
					std::vector<size_t> recog_neigs = { front_spaces[j] };
					//找最大的两个
					size_t num_neis = 2;
					auto temp_prob = space_prob;
					if (temp_prob.size() > 0) {
						size_t sele_num = temp_prob.size() > 2 ? num_neis : temp_prob.size();
						for (size_t k = 0; k < sele_num; ++k) {
							auto max_pair = *std::max_element(temp_prob.begin(), temp_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
							recog_neigs.push_back(max_pair.first);
							temp_prob[max_pair.first] = 0.;
						}
						m_recog_neighs.emplace_back(recog_neigs);
					}
				}
				for (size_t j = 0; j < exploit_num; ++j) {
					size_t space_inx;
					if (j < exploit_num / 2) {
						space_inx = std::distance(front_fre.begin(), std::min_element(front_fre.begin(), front_fre.end()));
						front_fre[space_inx]++;
						space_inx = std::distance(space_density.begin(), std::min_element(space_density.begin(), space_density.end()));
						space_density[space_inx] = (space_density[space_inx] * box_span[space_inx] + 1) / box_span[space_inx];
					}
					else {
						space_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
					}
					// method 1
					/*auto link_neis = getMO_HLC().getSubspaceInfo(front_spaces[space_inx]).m_link_subspaces;
					link_neis.push_back(front_spaces[space_inx]);*/
					//method 2
					//auto space_prob=spaceLinkProbability(front_spaces[space_inx]);
					auto space_prob = all_front_nei_prob[space_inx];
					size_t num_neis = 2;
					std::vector<size_t> neis;
					neis.push_back(front_spaces[space_inx]);
					for (size_t k = 0; k < num_neis; ++k) {
						Real rand_num = rnd->uniform.next();
						Real temp_prob = 0.;
						for (auto kk : space_prob) {
							if (rand_num > temp_prob && rand_num < temp_prob + kk.second) {
								if (std::find(neis.begin(), neis.end(), kk.first) == neis.end()) {
									neis.push_back(kk.first);
								}
								break;
							}
							else {
								temp_prob += kk.second;
							}
						}
					}
					auto off = sampleInFrontNeighSpace(neis, j, bound, 1, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
				
				operator_ratio[1] = 1.;
			}
			else if (type[i] == 3) {
				//k近邻交互
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					auto neighs = getPop()[i].nearestNeighbour(j,pro,m_neighs);
					std::vector<size_t> ind_inx;
					ind_inx.push_back(j);
					for (auto& nei : neighs) {
						ind_inx.push_back(nei.second);
					}
					SPMOEA_pop temp_pop(0, pro);
					for (size_t k = 0; k < ind_inx.size(); ++k) {
						temp_pop.append(getPop()[i][ind_inx[k]]);
					}
					auto off = sampleByDE(temp_pop, bound, 1, 2, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else if (type[i] == 4) {//基于个体爬山
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					auto off = samplePush(getPop()[i][j],1,pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else if (type[i] == 5) {//流形学习，umap降维

			}

			getOperatorRatio().emplace_back(operator_ratio);

			auto& interactive_sols = getInteractiveSols();
			for (size_t k = 0; k < interactive_sols.size(); ++k) {
				auto& ind = interactive_sols[k].back();
				getPop()[i].getOffspring()[k].variable() = ind->variable();
				getPop()[i].getOffspring()[k].objective() = ind->objective();
				getPop()[i].getOffspring()[k].setTimeEvaluate(0);
			}
		}
		updateEE(0, num_exploit);
	}

	bool SPMOEA1_17::spaceSubdivision(Problem *pro, Random *rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		
		auto front_spaces = getFrontSpace();
		//细分子空间
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			//size_t min_num = m_divide_granularity;
			size_t min_num = getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity;
			// min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
			//min_num = 20;
			if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
				//先得到子空间历史解用于更新细分的子空间的信息
				int dim = findSplitDim(sp, pro);
				auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
				Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
				splitSpace(sp, 2, dim, pos, false, pro, rnd);
			}
			else if (getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity<m_divide_granularity) {
				getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity++;
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

	void SPMOEA1_17::updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd) {
		auto& neigh = getMO_HLC().getSubspaceInfo(inx).m_sub_neighbors;
		std::vector<size_t> link_neighs;
		for (auto& i : neigh) {
			if (std::find(front_spaces.begin(), front_spaces.end(), i) != front_spaces.end()) {
				Real link_probability = ifContinuousSpace(inx, i);
				Real rand = rnd->uniform.next();
				if (rand <= link_probability) {
					link_neighs.push_back(i);
				}
			}
		}
		getMO_HLC().getSubspaceInfo(inx).m_link_subspaces = link_neighs;
	}

	void SPMOEA1_17::clusterSubspace() {
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

	std::vector<std::vector<size_t>> SPMOEA1_17::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_17::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA1_17::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_17::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem *pro, Random *rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA1_17::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_17::findSplitDim(int inx, Problem *pro) {
		////根据在哪一维上平分两边的子空间前沿解的差额最大来确定划分的维度
		//auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		//auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		//std::vector<Real> front_ind_ratio;
		//for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
		//	Real middle_v = (space_bound[j].first + space_bound[j].second) / 2;
		//	size_t upper_count = 0, lower_count = 0;
		//	for (size_t k = 0; k < front_sols.size(); ++k) {
		//		if (front_sols[k]->variable()[j] >= middle_v) {
		//			upper_count++;
		//		}
		//		else {
		//			lower_count++;
		//		}
		//	}
		//	size_t max_v = upper_count > lower_count ? upper_count : lower_count;
		//	size_t min_v = upper_count < lower_count ? upper_count : lower_count;
		//	if (max_v == 0) {
		//		max_v += 1;
		//	}
		//	front_ind_ratio.push_back((Real)min_v / max_v);
		//}
		//int dim = std::distance(front_ind_ratio.begin(), std::min_element(front_ind_ratio.begin(), front_ind_ratio.end()));

		//使用跨度确定分割的维度
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		int dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));
		
		////当维度跨度相差不大时，随机选择一维
		//Real max_span = *std::max_element(dim_span.begin(), dim_span.end());
		//std::vector<size_t> range_span_inx;
		//for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
		//	if (dim_span[j] / max_span >= 0.9) {
		//		range_span_inx.push_back(j);
		//	}
		//}
		//size_t select_inx = (size_t)std::floor(range_span_inx.size()*rnd->uniform.next());
		//int dim = range_span_inx[select_inx];
		return dim;
	}

	void SPMOEA1_17::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_17::recordMetrics(Problem *pro, Algorithm *alg) {
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

	void SPMOEA1_17::initiObjSpace(Problem *pro) {

	}
}