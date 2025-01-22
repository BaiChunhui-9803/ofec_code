#include "spmoea4_4.h"
#include "../../../../../utility/clustering/dbscan.h"
#include "../../../../../utility/clustering/kmeans.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA4_4::initialize_() {
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
		m_explore_rate = v.get<Real>("explore rate");
		m_switch_period = m_num_push + m_num_extend;
		//m_neighs = 5;
	}

	void SPMOEA4_4::run_() {
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

	void SPMOEA4_4::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA4_4::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(getPop().size());
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j].phenotype());
				}
			}
			m_off_solution.clear();
			m_off_solution.resize(getPop().size());
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_off_solution[i].push_back(&getPop()[i].getOffspring()[j].phenotype());
				}
			}
			m_his_solution.clear();
			auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_his_solution.push_back(&his_sols[i]->phenotype());
			}

			/*for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				auto& subspace_front_sols = getMO_HLC().getSubspaceInfo(i).m_subspace_front_sol;
				for (size_t j = 0; j < subspace_front_sols.size(); ++j) {
					m_solution.back().push_back(&subspace_front_sols[j]->phenotype());
				}
			}*/
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA4_4::initiVarSpace(Problem* pro) {
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

	void SPMOEA4_4::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
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

		size_t explore_pop_num = (size_t)std::floor(pop_num * m_explore_rate);
		size_t main_pop_num = pop_num - explore_pop_num;

		std::vector<size_t> inx1;
		for (size_t i = 0; i < main_pop_num; ++i) {
			inx1.push_back(i);
		}
		m_pop_index.emplace_back(inx1);
		std::vector<size_t> inx2;
		for (size_t i = main_pop_num; i < main_pop_num + explore_pop_num; ++i) {
			inx2.push_back(i);
		}
		m_pop_index.emplace_back(inx2);
		m_reset_box = bound;

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
		updateVarSpaceRank(pro, rnd);
		updateEE(pop_num, 0);
		SPMOEA::recordMetrics(pro, alg);

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

		m_operator_results;
		for (size_t i = 0; i < 8; ++i) {
			std::vector<size_t> temp = { 0,0,0 };
			m_operator_results.emplace_back(temp);
		}

		//SPMOEA::record();
	}

	int SPMOEA4_4::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
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
					Real rand = rnd->uniform.next();
					if (rand < m_explore_rate) {
						interactive_type.push_back(4);
					}
					else {
						interactive_type.push_back(7);
					}
				}
				else {
					interactive_type.push_back(8);
				}
				//interactive_type.push_back(6);
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
			for (size_t j = 0; j < getPop()[i].size(); j++) {
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
		//根据子代成功率自适应划分
		//adaptiveSplit(pro);
		auto pre_num_spaces = getMO_HLC().numSubspace();
		adaptiveSplitFrontSpace(pro, rnd);

		Real total_volume = getVarSpaceVolume();

		//更新子空间连续性
		auto cur_num_spaces = getMO_HLC().numSubspace();
		if (cur_num_spaces > pre_num_spaces) {
			updateFrontSpace();
			//更新子空间邻域
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
			}
			//更新种群分布和子代分布

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

		//testCoverage(pro);

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

		if (m_divide_iteration % m_switch_period < m_num_push) {
			//主种群
			SPMOEA_pop temp_pop1(0, pro);
			for (size_t j = 0; j < m_pop_index[0].size(); ++j) {
				temp_pop1.append(getPop()[0][m_pop_index[0][j]]);
			}
			for (size_t j = 0; j < m_pop_index[0].size(); ++j) {
				temp_pop1.append(getPop()[0].getOffspring()[m_pop_index[0][j]]);
			}
			SPMOEA::NDSort(temp_pop1);
			//探索种群
			SPMOEA_pop temp_pop2(0, pro);
			for (size_t j = 0; j < m_pop_index[1].size(); ++j) {
				temp_pop2.append(getPop()[0][m_pop_index[1][j]]);
			}
			for (size_t j = 0; j < m_pop_index[1].size(); ++j) {
				temp_pop2.append(getPop()[0].getOffspring()[m_pop_index[1][j]]);
			}
			SPMOEA::NDSort(temp_pop2);

			std::vector<size_t> keep_inx;
			for (size_t i = 0; i < temp_pop2.size(); ++i) {
				auto obj1 = temp_pop2[i].objective();
				bool flag = true;
				for (size_t j = 0; j < temp_pop1.size(); ++j) {
					if (temp_pop1[j].fitness() == 0) {
						auto obj2 = temp_pop1[j].objective();
						auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
						if (ship == Dominance::kDominated) {
							flag = false;
							break;
						}
					}
				}
				if (flag) {
					keep_inx.push_back(i);
				}
			}
			//更新主种群
			/*for (size_t j = 0; j < keep_inx.size(); ++j) {
				temp_pop1.append(temp_pop2[keep_inx[j]]);
			}
			SPMOEA::NDSort(temp_pop1);*/
			std::vector<std::vector<Real>> all_data;
			std::vector<size_t> ranks;
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				all_data.emplace_back(temp_pop1[j].objective());
				ranks.push_back(temp_pop1[j].fitness());
			}
			size_t sele_num = m_pop_index[0].size()-keep_inx.size();
			auto sele_inx = selectMaxMinFromFront(all_data, ranks, sele_num);
			for (size_t j = 0; j < m_pop_index[0].size()- keep_inx.size(); ++j) {
				getPop()[0][m_pop_index[0][j]] = temp_pop1[sele_inx[j]];
			}
			for (size_t j = 0; j < keep_inx.size(); ++j) {
				getPop()[0][sele_num+j] = temp_pop2[keep_inx[j]];
			}

			//更新探索种群，子空间内排序，与历史解一起排
			if (keep_inx.size() >= 0.5*m_pop_index[1].size() || m_explore_fail_times >= 10) {
				m_explore_reset_flag = true;
				m_explore_fail_times = 0;
				for (size_t j = 0; j < m_pop_index[1].size(); ++j) {
					getPop()[0][m_pop_index[1][j]] = getPop()[0].getOffspring()[m_pop_index[1][j]];
				}
			}
			else {
				if (m_explore_reset_flag) {
					for (size_t j = 0; j < m_pop_index[1].size(); ++j) {
						getPop()[0][m_pop_index[1][j]] = getPop()[0].getOffspring()[m_pop_index[1][j]];
					}
					m_explore_reset_flag = false;
				}
				else {
					std::vector<std::vector<Real>> all_data2;
					std::vector<size_t> ranks2;
					for (size_t j = 0; j < temp_pop2.size(); ++j) {
						all_data2.emplace_back(temp_pop2[j].objective());
						ranks2.push_back(temp_pop2[j].fitness());
					}
					auto sele_inx2 = selectMaxMinFromFront(all_data2, ranks2, m_pop_index[1].size());
					for (size_t j = 0; j < m_pop_index[1].size(); ++j) {
						getPop()[0][m_pop_index[1][j]] = temp_pop2[sele_inx2[j]];
					}
					size_t better_off_num = 0;
					size_t first_num = 0;
					for (size_t j = 0; j < temp_pop2.size(); ++j) {
						if (temp_pop2[j].fitness() == 0) {
							first_num++;
							if (j >= m_pop_index[1].size()) {
								better_off_num++;
							}
						}
					}
					if (keep_inx.size()<1 && (Real)better_off_num/first_num<0.3) {
						m_explore_fail_times++;
					}
					else{
						m_explore_fail_times = 0;
					}

					std::cout << "the number of front solutions: " << keep_inx.size() << std::endl;
					std::cout << "the  number   of  better offs: " << better_off_num << std::endl;
					//std::cout << "the  number   of  first offs: " << first_num << std::endl;
				}
			}
		}

		m_subspace_pop_index.clear();

		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA4_4::assignPop(Population<Solution<>>& pop, Problem* pro, Algorithm* alg, Random* rnd) {
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
			Real v = getMO_HLC().subspaceTree().getBoxVolume(behind_spaces[i]);
			box_span.push_back(v);
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
			std::vector<Real> temp = { space_density[i],space_rank[i] };
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//第一层子空间
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(behind_spaces[i]);
			}
		}
		//随机选择一个子空间
		if (first_layer.size() > 0 && pop.size() > 0) {
			size_t space_idx = (size_t)std::floor(first_layer.size() * rnd->uniform.next());
			auto& his_sols = getMO_HLC().getSubspaceInfo(first_layer[space_idx]).m_history_inds;
			size_t count = 0;
			if (his_sols.size() > pop.size()) {
				std::vector<size_t> sele_inx;
				while (sele_inx.size() < pop.size()) {
					size_t idx = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					if (sele_inx.empty()) {
						sele_inx.push_back(idx);
					}
					else if (std::find(sele_inx.begin(), sele_inx.end(), idx) == sele_inx.end()) {
						sele_inx.push_back(idx);
					}
				}
				for (size_t i = 0; i < sele_inx.size(); ++i) {
					pop[i].variable() = his_sols[sele_inx[i]]->variable();
					pop[i].objective() = his_sols[sele_inx[i]]->objective();
				}
			}
			else {
				while (count < pop.size()) {
					size_t inx = (size_t)std::floor(behind_spaces.size() * rnd->uniform.next());
					if (getMO_HLC().getSubspaceInfo(behind_spaces[inx]).m_history_inds.size() > 0) {
						auto& his_sol = getMO_HLC().getSubspaceInfo(behind_spaces[inx]).m_history_inds;
						size_t idx = (size_t)std::floor(his_sol.size() * rnd->uniform.next());
						pop[count].variable() = his_sol[idx]->variable();
						pop[count].objective() = his_sol[idx]->objective();
						count++;
					}
				}


			}
		}

		m_immigrant_times += pop.size();
	}

	void SPMOEA4_4::testCoverage(Problem* pro) {
		//检测子空间是否包含不同的PS片段
		m_real_front_subspaces.clear();
		m_multi_segment_subspaces.clear();
		m_match_subspace.clear();
		m_error_subspace.clear();
		m_loss_subspace.clear();
		m_match_multi_seg_subspace.clear();
		auto front_sp = getFrontSpace();
		//动态分割下使用下面的代码，真实PS所在的子空间
		//std::vector<size_t> real_front_subspaces;
		////具有多段的子空间数目
		//std::vector<size_t> multi_segment_subspaces;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			auto& box = getMO_HLC().subspaceTree().getBox(i);
			//沿着第一维构造解
			size_t num = 200;
			std::vector<size_t> flag(num - 1, 1);
			Real delta = (box[0].second - box[0].first) / num;
			for (size_t j = 0; j < num - 1; ++j) {
				std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
				auto sol = CAST_CONOP(pro)->createVar(s);
				//判断sol是否在box内
				for (size_t k = 0; k < box.size(); ++k) {
					if (sol[k] > box[k].second || sol[k] < box[k].first) {
						flag[j] = 0;
						break;
					}
				}
			}
			//统计段数
			size_t segment_time = 0;
			for (size_t j = 1; j < flag.size(); ++j) {
				if (flag[j - 1] == 0 && flag[j] == 1) {
					segment_time++;
				}
			}
			if (flag[0] == 1) {
				segment_time++;
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

		////测试计算邻域识别的准确性
		//if (m_divide_iteration % m_switch_period >= m_num_push) {
		//	std::cout << "************************************************************** " << std::endl;
		//	std::cout << "neighs match result: " << std::endl;
		//	int total_match_num = 0;
		//	for (size_t i = 0; i < m_recog_neighs.size(); ++i) {
		//		int match_num = 0;
		//		auto center_space = m_recog_neighs[i][0];
		//		if (std::find(m_real_front_subspaces.begin(), m_real_front_subspaces.end(), center_space) != m_real_front_subspaces.end()) {
		//			//找哪些位置匹配
		//			std::vector<size_t> match_inx;
		//			for (size_t j = 0; j < m_real_all_front_subspaces.size(); ++j) {
		//				if (m_real_all_front_subspaces[j] == center_space) {
		//					match_inx.push_back(j);
		//				}
		//			}
		//			for (size_t j = 0; j < match_inx.size(); ++j) {
		//				if (match_inx[j] == 0) {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] + 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1) {
		//							match_num++;
		//						}
		//					}
		//				}
		//				else if (match_inx[j] == m_real_all_front_subspaces.size() - 1) {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1) {
		//							match_num++;
		//						}
		//					}
		//				}
		//				else {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
		//					size_t pos2 = m_real_all_front_subspaces[match_inx[j] + 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1 || m_recog_neighs[i][k] == pos2) {
		//							match_num++;
		//						}
		//					}
		//				}
		//			}
		//			total_match_num += match_num;
		//		}
		//		else {
		//			match_num = -1;
		//		}
		//		std::cout << match_num << " ";
		//	}
		//	std::cout << std::endl;
		//	//输出总的匹配率
		//	std::cout << "neighs match result: " << std::setprecision(2) << (Real)total_match_num / (2. * m_real_front_subspaces.size() - 2) << std::endl;
		//}

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
			auto& his_sols = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_history_inds;
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
			std::cout << std::setprecision(2) << "(" << span[min_inx] << "," << span[max_inx] << ")";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "the   number   of   immigrant  is: " << m_immigrant_times << std::endl;
		std::cout << std::endl;

		std::cout << std::endl;
	}

	//子空间的前沿子空间邻域交互形成一个解当做子代
	std::vector<std::vector<Real>> SPMOEA4_4::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		//子连通空间内采样
		int method = 1;
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		size_t con_num = CAST_CONOP(pro)->numberConstraints();
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			auto sol_num = getMO_HLC().getSubspaceInfo(front_link_spaces[i]).m_history_inds.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num / volume);
		}
		auto& sol = getPop()[0][ind_inx].variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//选出一个子空间
			size_t sele_inx = 0;
			if (method == 1) {//随机选择
				sele_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//根据子空间内的前沿数选择
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				//front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[inx]);
			}
			sele_inx = front_link_spaces[0];
			//找出其前沿邻域子空间
			std::vector<size_t> nei1;
			/*nei1.push_back(front_link_spaces[sele_inx]);
			auto neigh = getMO_HLC().getSubspaceInfo(front_link_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_link_spaces.begin(), front_link_spaces.end(), jj) != front_link_spaces.end()) {
					nei1.push_back(jj);
				}
			}*/
			for (auto jj : front_link_spaces) {
				nei1.push_back(jj);
			}
			//将邻域前沿个体构成新种群
			size_t pop_size = 0;
			size_t first_space_num = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				pop_size += front_sol.size();
				if (j == 1) {
					first_space_num = front_sol.size();
				}
			}
			PopDE<> temp_pop(pop_size, pro);
			size_t count = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				for (size_t k = 0; k < front_sol.size(); ++k) {
					temp_pop[count].variable() = front_sol[k]->variable();
					temp_pop[count].objective() = front_sol[k]->objective();
					count++;
				}
			}
			std::vector<size_t> candidate_inx;
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(first_space_num * rnd->uniform.next());
				candidate_inx.push_back(inx);
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				temp_pop[candidate_inx[0]].donor().variable() = off1[0];
				temp_pop[candidate_inx[0]].donor().objective() = getInteractiveSols().back().back()->objective();
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//子空间内随机交互
				candidate_inx.push_back((size_t)std::floor(first_space_num * rnd->uniform.next()));
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[1]);
				candidate_inx.push_back(inx);
				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);

				temp_pop[candidate_inx[0]].donor().evaluate(pro, this);
				Solution<> ind4(temp_pop[candidate_inx[0]].donor());
				/*ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();*/
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].donor().variable().vect());
			}
			m_total_extend_times++;
			//父代与子代的支配关系
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].donor().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_extend_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_extend_results[0]++;
			}
			else {
				m_extend_results[2]++;
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_4::sampleLinearInFrontSpace(size_t space, size_t ind_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		//子连通空间内采样
		auto& sol = getPop()[0][ind_inx].variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//将邻域前沿个体构成新种群
			size_t pop_size = 0;
			auto& front_sol = getMO_HLC().getSubspaceInfo(space).m_front_sol_in_subspace;
			pop_size += front_sol.size();
			PopDE<> temp_pop(pop_size, pro);
			for (size_t k = 0; k < front_sol.size(); ++k) {
				temp_pop[k].variable() = front_sol[k]->variable();
				temp_pop[k].objective() = front_sol[k]->objective();
			}
			std::vector<size_t> candidate_inx;
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(pop_size * rnd->uniform.next());
				candidate_inx.push_back(inx);
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				temp_pop[candidate_inx[0]].donor().variable() = off1[0];
				temp_pop[candidate_inx[0]].donor().objective() = getInteractiveSols().back().back()->objective();
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//子空间内向量采样
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[0]);
				candidate_inx.push_back(inx);

				auto sol1 = temp_pop[candidate_inx[0]].variable().vect();
				auto sol2 = temp_pop[candidate_inx[1]].variable().vect();
				Real rand = rnd->uniform.next();

				std::vector<Real> new_sol;
				for (size_t k = 0; k < sol1.size(); ++k) {
					new_sol.push_back(sol1[k] + rand * (sol2[k] - sol1[k]));
				}
				temp_pop[candidate_inx[0]].donor().variable().vect() = new_sol;
				temp_pop[candidate_inx[0]].donor().evaluate(pro, this);
				Solution<> ind4(temp_pop[candidate_inx[0]].donor());

				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);

				out_off.emplace_back(temp_pop[candidate_inx[0]].donor().variable().vect());
			}
			m_total_extend_times++;
			//父代与子代的支配关系
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].donor().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_extend_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_extend_results[0]++;
			}
			else {
				m_extend_results[2]++;
			}
		}
		return out_off;
	}

	//基于子空间个体采样，向个体所在子空间前沿解推进，提升个体
	std::vector<std::vector<Real>> SPMOEA4_4::samplePush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<size_t> equal_space;
		auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == getMO_HLC().getSubspaceInfo(space).m_best_rank && getMO_HLC().getSubspaceInfo(space).m_best_rank != INT16_MAX) {
				equal_space.push_back(ii);
			}
		}
		std::vector<size_t> better_space;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < getMO_HLC().getSubspaceInfo(space).m_best_rank) {
				better_space.push_back(ii);
			}
		}
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(sol);
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable().vect() = off[0];
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//根据所在子空间位置概率选择推进方向
				auto base_sol = sol.variable().vect();
				//先看该解是不是子空间前沿解
				size_t sol_inx = 0;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto temp_sol = his_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						sol_inx = j;
						break;
					}
				}
				//子空间历史解分两类：前沿和非前沿
				std::vector<size_t> subspace_front_inx;
				std::vector<size_t> subspace_behind_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					if (his_sols[j]->type() == 0) {//type==0, is front sol
						subspace_front_inx.push_back(j);
					}
					else {
						subspace_behind_inx.push_back(j);
					}
				}
				Real rand = rnd->uniform.next();
				temp_pop.append(*his_sols[sol_inx]);
				candidate_inx = { 0,1,2 };
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					if (rand < push_prob) {
						if (better_space.empty()) {
							//局部最优子空间，内部交互
							if (subspace_behind_inx.empty()) {//自由交互
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
								m_push_behind_empty++;
							}
							else {
								//子空间内前后排交互
								size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);

								////子空间内找支配对
								//size_t a = m_evaluations;
								//size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								//std::vector<size_t> dominate_inx;
								//auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								//for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
								//	auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
								//	Dominance ship = objectiveCompare(p2_obj,p1_obj,pro->optimizeMode());
								//	if (ship == Dominance::kDominant) {
								//		dominate_inx.push_back(j);
								//	}
								//}
								//size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								//temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);
								//temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);

							}
						}
						else {//选择邻域一个更好子空间的前沿解
							size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
							//邻域子空间找支配对

						}
					}
					else {
						if (equal_space.empty()) {
							if (subspace_front_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[inx1]]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[subspace_front_inx[inx2]]);
							}
							/*size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							temp_pop.append(*his_sols[inx1]);
							size_t inx2 = 0;
							do {
								inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							} while (inx2 == inx1);
							temp_pop.append(*his_sols[inx2]);*/
						}
						else {
							size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
						}
					}

				}
				else {//选择子空间中支配其的前沿解作差分
					Real external_prob = 0.5;
					if (rand < 1 - push_prob) {//前推，与前排解交互
						size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);
						size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);

						////子空间内找支配对
						//size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//std::vector<size_t> dominate_inx;
						//auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
						//for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
						//	auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
						//	Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
						//	if (ship == Dominance::kDominant) {
						//		dominate_inx.push_back(j);
						//	}
						//}
						//size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
						//temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

						//temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
					}
					else {//与邻域更好的子空间交互
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[inx1]);
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(*his_sols[inx2]);
						//if (better_space.empty()) {//局部最优子空间，内部交互
						//	candidate_inx = { 0,1,2 };
						//	/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

						//	size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);*/

						//	//子空间内找支配对
						//	size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//	std::vector<size_t> dominate_inx;
						//	auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
						//	for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
						//		auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
						//		Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
						//		if (ship == Dominance::kDominant) {
						//			dominate_inx.push_back(j);
						//		}
						//	}
						//	size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

						//	temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
						//}
						//else {//选择邻域一个更好子空间的前沿解
						//	candidate_inx = { 0,1,2 };
						//	size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
						//	size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
						//	temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

						//	size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
						//}
					}
				}
				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
			m_total_push_times++;
			//父代与子代的支配关系
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].trial().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_push_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_push_results[0]++;
			}
			else {
				m_push_results[2]++;
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_4::sampleByIndPos(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<size_t> equal_space;
		auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
		size_t space_rank = getMO_HLC().getSubspaceInfo(space).m_best_rank;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == space_rank && space_rank != INT16_MAX) {
				equal_space.push_back(ii);
			}
		}
		std::vector<size_t> better_space;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < space_rank) {
				better_space.push_back(ii);
			}
		}
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			size_t flag = 0;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(sol);
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable() = sol.variable();
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//根据所在子空间位置概率选择推进方向
				auto base_sol = sol.variable().vect();
				//先看该解是不是子空间前沿解
				size_t sol_inx = 0;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto temp_sol = his_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						sol_inx = j;
						break;
					}
				}
				//子空间历史解分两类：前沿和非前沿
				std::vector<size_t> subspace_front_inx;
				std::vector<size_t> subspace_behind_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					if (his_sols[j]->type() == 0) {//type==0, is front sol
						subspace_front_inx.push_back(j);
					}
					else {
						subspace_behind_inx.push_back(j);
					}
				}
				Real rand = rnd->uniform.next();
				temp_pop.append(*his_sols[sol_inx]);
				candidate_inx = { 0,1,2 };
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					Real rand2 = rnd->uniform.next();
					if (rand < push_prob) {//前推
						if (rand2 < 0.5) {//内部前推
							flag = 1;
							if (subspace_behind_inx.empty()) {//自由交互
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
								m_push_behind_empty++;
							}
							else {
								//子空间内前后排交互
								/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);
								*/
								//子空间内找支配对
								size_t a = m_evaluations;
								size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								std::vector<size_t> dominate_inx;
								auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
									auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
									Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
									if (ship == Dominance::kDominant) {
										dominate_inx.push_back(j);
									}
								}
								size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
							}
						}
						else {//外部前推
							flag = 2;
							if (better_space.empty()) {
								//局部最优子空间，内部交互
								if (subspace_behind_inx.empty()) {//自由交互
									size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[inx1]);
									size_t inx2 = 0;
									do {
										inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									} while (inx2 == inx1);
									temp_pop.append(*his_sols[inx2]);
									m_push_behind_empty++;
								}
								else {
									//子空间内前后排交互
									/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

									size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);
									*/
									//子空间内找支配对
									size_t a = m_evaluations;
									size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
									std::vector<size_t> dominate_inx;
									auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
									for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
										auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
										Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
										if (ship == Dominance::kDominant) {
											dominate_inx.push_back(j);
										}
									}
									size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

									temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
								}
							}
							else {//选择邻域一个更好子空间的前沿解
								size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_history_inds.size() * rnd->uniform.next());
								temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_history_inds[sele_ind_inx]);

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
								//邻域子空间找支配对

							}
						}
					}
					else {//扩展
						if (rand2 < 0.5) {//内部扩展
							flag = 3;
							if (subspace_front_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[inx1]]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[subspace_front_inx[inx2]]);
							}
						}
						else {//外部扩展
							flag = 4;
							if (equal_space.empty()) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
							}
							else {
								size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
							}
						}
					}
				}
				else {//选择子空间中支配其的前沿解作差分
					Real rand2 = rnd->uniform.next();
					if (rand < 1 - push_prob) {//前推，与前排解交互
						if (rand2 < 0.5) {//内部前推
							flag = 5;
							size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
							std::vector<size_t> dominate_inx;
							auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
							for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
								auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
								Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
								if (ship == Dominance::kDominant) {
									dominate_inx.push_back(j);
								}
							}
							size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
							temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

							temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
						}
						else {//外部前推
							flag = 6;
							if (better_space.empty()) {//局部最优子空间，内部交互
								/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);*/

								//子空间内找支配对
								size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								std::vector<size_t> dominate_inx;
								auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
									auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
									Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
									if (ship == Dominance::kDominant) {
										dominate_inx.push_back(j);
									}
								}
								size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
							}
							else {//选择邻域一个更好子空间的前沿解
								size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
							}
						}
					}
					else {
						if (rand2 < 0.5) {//内部扩展
							flag = 7;
							if (subspace_behind_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[inx1]]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[subspace_behind_inx[inx2]]);
							}
						}
						else {//外部扩展
							flag = 8;
							if (equal_space.empty()) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[inx1]);
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(*his_sols[inx2]);
							}
							else {
								size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_history_inds.size() * rnd->uniform.next());
								temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_history_inds[sele_ind_inx]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);
							}
						}
					}
				}
				temp_pop[candidate_inx[0]].mutate(temp_pop.scalingFactor(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
			if (flag != 0) {
				m_ind_operator_times[flag - 1]++;
				//父代与子代的支配关系
				auto& pp = temp_pop[candidate_inx[0]].objective();
				auto& off = temp_pop[candidate_inx[0]].trial().objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					m_operator_results[flag - 1][1]++;
				}
				else if (dominance_ship == Dominance::kDominated) {
					m_operator_results[flag - 1][0]++;
				}
				else {
					m_operator_results[flag - 1][2]++;
				}
			}

		}
		return out_off;
	}

	//个体层面的选择，非支配时随机选择，子代更差时，根据父代是否为前沿决定子代保留的概率
	void SPMOEA4_4::SolutionSelection(Problem* pro, Random* rnd) {
		//后续解的选择需要考虑离已选解的分布
		//综合当前种群、子代种群和历史非支配解在linkspace的分布进行选择
		//子空间先随机选一个，然后依次选择离linkspace子空间已选解最远的一个候选解
		auto& his_front_sols = getHisFrontSols();
		//根据有解分布的子空间是否为前沿子空间来分配，先在各个子空间随机选择一个当前靠前的解，然后选择非前沿子空间中较好的解，
		// 最后轮流依次遍历前沿子空间的linkspace,选择与已选解最远的一个

		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			Population<Solution<>> residual_pop;
			std::vector<size_t> sele_inx;
			std::vector<size_t> residual_inx;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				bool flag1 = false;//父代个体是否为前沿个体
				auto& sol1 = getPop()[i][j].variable().vect();
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					auto sol2 = his_front_sols[k]->variable().vect();
					if (ifSame(sol1, sol2)) {
						flag1 = true;
						break;
					}
				}
				bool flag2 = false;//子代个体是否为前沿个体
				auto& sol2 = getPop()[i][j].variable().vect();
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					auto sol3 = his_front_sols[k]->variable().vect();
					if (ifSame(sol2, sol3)) {
						flag2 = true;
						break;
					}
				}

				if (dominance_ship == Dominance::kDominant) {//父代支配子代
					//根据个体是否为前沿解决定替换概率
					//Real rand = rnd->uniform.next();
					/*if (flag1) {
						if (rand < 0.1) {
							sele_inx.push_back(j);
						}
						else {
							sele_inx.push_back(getPop()[i].size() + j);
						}
					}
					else {
						if (rand < 0.2) {
							sele_inx.push_back(j);
						}
						else {
							sele_inx.push_back(getPop()[i].size() + j);
						}
					}*/
					sele_inx.push_back(getPop()[i].size() + j);
					m_total_results[0]++;
				}
				else if (dominance_ship == Dominance::kDominated) {//子代支配父代
					/*if (flag2) {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					else {
						sele_inx.push_back(j);
					}*/
					sele_inx.push_back(j);
					m_total_results[1]++;
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						sele_inx.push_back(j);
					}
					else {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					//auto pp_sol = getPop()[i][j].variable().vect();
					//size_t p_space = getMO_HLC().subspaceTree().getRegionIdx(pp_sol);
					//auto off_sol = getPop()[i].getOffspring()[j].variable().vect();
					//size_t o_space = getMO_HLC().subspaceTree().getRegionIdx(off_sol);
					//Real p_span = 0.;
					//Real o_span = 0.;
					//auto p_box = getMO_HLC().subspaceTree().getBox(p_space);
					//for (size_t j = 0; j < p_box.size(); ++j) {
					//	p_span += (p_box[j].second - p_box[j].first);
					//}
					//p_span /= p_box.size();
					//auto o_box = getMO_HLC().subspaceTree().getBox(o_space);
					//for (size_t j = 0; j < o_box.size(); ++j) {
					//	o_span += (o_box[j].second - o_box[j].first);
					//}
					//o_span /= o_box.size();

					//if (flag1 && flag2) {
					//	//比较前沿解的个数
					//	Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_front_sol_in_subspace.size() / p_span;
					//	Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_front_sol_in_subspace.size() / o_span;
					//	if (density1 < density2) {
					//		sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else {
					//		sele_inx.push_back(j);
					//	}
					//}
					//else if (flag1 && (!flag2)) {
					//	sele_inx.push_back(j);
					//}
					//else if (flag2 && (!flag1)) {
					//	sele_inx.push_back(getPop()[i].size() + j);
					//}
					//else {
					//	/*Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / p_span;
					//	Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / o_span;
					//	if (density1 < density2) {
					//		sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else {
					//		sele_inx.push_back(j);
					//	}*/
					//	size_t rank1 = getMO_HLC().getSubspaceInfo(p_space).m_best_rank;
					//	size_t rank2 = getMO_HLC().getSubspaceInfo(o_space).m_best_rank;
					//	if (rank1 < rank2) {
					//		Real rand = rnd->uniform.next();
					//		if (rand < 0.2) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//		//sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else if (rank1 > rank2) {
					//		Real rand = rnd->uniform.next();
					//		if (rand > 0.2) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//		//sele_inx.push_back(j);
					//	}
					//	else {
					//		Real rand = rnd->uniform.next();
					//		if (rand > 0.5) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//	}
					//}
					m_total_results[2]++;
				}
			}
			//使用选出的新解更新个体演化轨迹
			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}


	void SPMOEA4_4::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA4_4::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA4_4::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
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

	void SPMOEA4_4::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		//根据子空间的排序、前沿解的跨度
		auto front_spaces = getFrontSpace();

		assign_pop_resource.push_back(max_pop_size);
		//updatePopResource(assign_pop_resource);
	}

	void SPMOEA4_4::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
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
			//m_recog_neighs.clear();
			//m_all_front_nei_prob.clear();
			//std::vector<std::map<size_t, Real>> all_front_nei_prob;
			//for (size_t j = 0; j < front_spaces.size(); ++j) {
			//	auto space_prob = frontSpaceLinkProbability(front_spaces[j]);
			//	all_front_nei_prob.emplace_back(space_prob);
			//	/***************测试识别的邻域关系是否准确**************/
			//	std::vector<size_t> recog_neigs = { front_spaces[j] };
			//	//找最大的两个
			//	size_t num_neis = 2;
			//	auto temp_prob = space_prob;
			//	if (temp_prob.size() > 0) {
			//		size_t sele_num = temp_prob.size() > 2 ? num_neis : temp_prob.size();
			//		for (size_t k = 0; k < sele_num; ++k) {
			//			auto max_pair = *std::max_element(temp_prob.begin(), temp_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
			//			recog_neigs.push_back(max_pair.first);
			//			temp_prob[max_pair.first] = 0.;
			//		}
			//		m_recog_neighs.emplace_back(recog_neigs);
			//	}
			//}
			//m_all_front_nei_prob = all_front_nei_prob;
			if (type[i] == 0) {
				//随机初始化
				all_off = sampleRandom(getPop()[i], search_bound, exploit_num, pro, alg, rnd);
			}
			else if (type[i] == 1) {
				//auto all_off = sampleByGA(getPop()[i], exploit_num, pro, rnd);
				size_t kk = 2;
				//all_off = sampleByDE(getPop()[i], bound, exploit_num, kk, pro, alg, rnd);
				all_off = sampleDiverse(getPop()[i], exploit_num, pro, alg, rnd);
			}
			else if (type[i] == 4) {//基于个体爬山
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					//auto off = samplePush(getPop()[i][j], 1, 0.2, pro, alg, rnd);
					//auto off = sampleByIndPos(getPop()[i][j], 1,0.5, pro, alg, rnd);
					auto off = sampleByPush(getPop()[i][j], 1, 0.5, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else if (type[i] == 7) {//判断邻域，邻域交互
				//个体所在子空间，与哪个子空间为邻域
				std::map<size_t, std::vector<size_t>> pop2space;
				m_subspace_pop_index.clear();
				m_pop_index_clusters.clear();
				auto first_pop = m_pop_index[0];
				auto second_pop = m_pop_index[1];
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					if (std::find(first_pop.begin(), first_pop.end(), j) != first_pop.end()) {
						auto var = getPop()[i][j].variable().vect();
						auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
						if (pop2space[space_idx].empty()) {
							std::vector<size_t> temp_inx;
							pop2space.insert(std::make_pair(space_idx, temp_inx));
						}
						pop2space[space_idx].emplace_back(j);
					}
				}
				m_subspace_pop_index = pop2space;
				auto front_spaces = getFrontSpace();
				//子空间内演化，和当前个体交互还是历史前沿交互？
				std::vector<std::pair<size_t, size_t>> add_pairs;
				for (auto& sub : pop2space) {
					size_t idx = sub.first;
					PopDE<> temp_pop;
					std::vector<size_t> candidate_inx;
					int method = 1;
					Real rand = rnd->uniform.next();
					if (rand < 0.5) {
						method = 1;
					}
					else {
						method = 2;
					}
					std::vector<size_t> act_space = { idx };
					std::vector<size_t> neis;
					if (method == 1) {//有个体的子空间
						for (auto& sub_nei : pop2space) {
							if (sub_nei.first != idx) {
								neis.push_back(sub_nei.first);
							}
						}
						auto nei_space_prob = spaceLinkProbability(idx, neis);
						std::vector<size_t> pop_index;
						if (nei_space_prob.size() > 0) {
							//直接取最大值
							auto max_pair = *std::max_element(nei_space_prob.begin(), nei_space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
							//和邻域子空间内的个体交互
							act_space.push_back(max_pair.first);
							if (nei_space_prob.size() > 1) {
								nei_space_prob[max_pair.first] = 0.;
								auto second_pair = *std::max_element(nei_space_prob.begin(), nei_space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
								Real rand = rnd->uniform.next();
								if (rand < 0.5) {
									act_space.back() = second_pair.first;
								}
							}
							for (size_t j = 0; j < act_space.size(); ++j) {
								for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
									pop_index.push_back(pop2space[act_space[j]][k]);
								}
							}
							std::pair<size_t, size_t> neigh_pair;
							neigh_pair.first = idx;
							neigh_pair.first = max_pair.first;
							bool flag = false;
							for (size_t j = 0; j < add_pairs.size(); ++j) {
								if ((neigh_pair.first == add_pairs[j].first && neigh_pair.second == add_pairs[j].second) || (neigh_pair.first == add_pairs[j].second && neigh_pair.second == add_pairs[j].first)) {
									flag = true;
								}
							}
							if (!flag) {
								add_pairs.emplace_back(neigh_pair);
								m_pop_index_clusters.emplace_back(pop_index);
							}
						}
						else {
							for (size_t j = 0; j < act_space.size(); ++j) {
								for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
									pop_index.push_back(pop2space[act_space[j]][k]);
								}
							}
							m_pop_index_clusters.emplace_back(pop_index);
						}
						for (size_t j = 0; j < act_space.size(); ++j) {
							for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
								temp_pop.append(getPop()[i][pop2space[act_space[j]][k]]);
							}
						}
					}
					else {//找邻域
						auto neighs = getMO_HLC().getSubspaceInfo(idx).m_sub_neighbors;
						for (auto& sub_nei : neighs) {
							neis.push_back(sub_nei);
						}
						auto nei_space_prob = spaceLinkProbability(idx, neis);
						if (nei_space_prob.size() > 0) {
							//直接取最大值
							auto max_pair = *std::max_element(nei_space_prob.begin(), nei_space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
							//和邻域子空间内的个体交互
							act_space.push_back(max_pair.first);
							if (nei_space_prob.size() > 1) {
								nei_space_prob[max_pair.first] = 0.;
								auto second_pair = *std::max_element(nei_space_prob.begin(), nei_space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
								Real rand = rnd->uniform.next();
								if (rand < 0.5) {
									act_space.back() = second_pair.first;
								}
								/*if (rand < second_pair.second / (max_pair.second + second_pair.second)) {
									act_space.back() = second_pair.first;
								}*/
							}
							//邻域是否有个体
							if (pop2space.find(act_space.back())==pop2space.end()) {
								for (size_t j = 0; j < act_space.size() - 1; ++j) {
									for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
										temp_pop.append(getPop()[i][pop2space[act_space[j]][k]]);
									}
								}
							}
							else {
								for (size_t j = 0; j < act_space.size(); ++j) {
									for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
										temp_pop.append(getPop()[i][pop2space[act_space[j]][k]]);
									}
								}
							}
						}
						else {
							for (size_t j = 0; j < act_space.size(); ++j) {
								for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
									temp_pop.append(getPop()[i][pop2space[act_space[j]][k]]);
								}
							}
						}
						for (size_t j = 0; j < act_space.size(); ++j) {
							if (std::find(front_spaces.begin(), front_spaces.end(), act_space[j]) == front_spaces.end()) {
								auto& front_sols = getMO_HLC().getSubspaceInfo(act_space[j]).m_subspace_front_sol;
								for (size_t j = 0; j < front_sols.size(); ++j) {
									temp_pop.append(*front_sols[j]);
								}
							}
							else {
								auto& front_sols = getMO_HLC().getSubspaceInfo(act_space[j]).m_front_sol_in_subspace;
								for (size_t j = 0; j < front_sols.size(); ++j) {
									temp_pop.append(*front_sols[j]);
								}
							}
						}
					}

					//产生一定个数的解
					for (size_t j = 0; j < sub.second.size(); ++j) {
						candidate_inx.clear();
						//size_t sele_inx1 = (size_t)std::floor(temp_pop.size()*rnd->uniform.next());
						candidate_inx.push_back(j);

						size_t sele_inx2 = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
						auto& sol2 = temp_pop[sele_inx2].variable().vect();
						candidate_inx.push_back(sele_inx2);
						size_t sele_inx3 = 0;
						std::vector<Real> sol3;
						do {
							sele_inx3 = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
							sol3 = temp_pop[sele_inx3].variable().vect();
						} while (ifSame(sol2, sol3));
						candidate_inx.push_back(sele_inx3);

						temp_pop[j].mutate(temp_pop.scalingFactor(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
						temp_pop.recombine(j, rnd, pro);
						temp_pop[j].trial().evaluate(pro, this);

						Solution<> ind4(temp_pop[j]);
						ind4.variable() = temp_pop[j].trial().variable();
						ind4.objective() = temp_pop[j].trial().objective();
						std::vector<std::shared_ptr<Solution<>>> temp_pair;
						temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[j]));
						temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
						temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
						temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
						getInteractiveSols().emplace_back(temp_pair);

						getPop()[i].getOffspring()[sub.second[j]].variable() = ind4.variable();
						getPop()[i].getOffspring()[sub.second[j]].objective() = ind4.objective();
					}
				}

				//探索种群在选定位置进行演化
				if (m_explore_reset_flag) {
					//选择一个子空间初始化
					//size_t space = (size_t)std::floor(total_spaces*rnd->uniform.next());
					auto sol_info = selectSubspaceFromObj(pro, rnd);
					//选择最大的一对点或者随机选择一对
					auto max_pair = *std::max_element(sol_info.begin(), sol_info.end(), [](std::pair<Real, std::pair<size_t, size_t>> left, std::pair<Real, std::pair<size_t, size_t>> right) {return left.first < right.first; });
					
					//根据点对确定bound
					auto point_pair = max_pair.second;
					size_t num_var = pro->numberVariables();
					auto sol1 = getPop()[i][first_pop[point_pair.first]].variable().vect();
					auto sol2 = getPop()[i][first_pop[point_pair.second]].variable().vect();
					
					/*auto& his_front_sols = getHisFrontSols();
					auto sol1 = his_front_sols[point_pair.first]->variable().vect();
					auto sol2 = his_front_sols[point_pair.second]->variable().vect();
					*/
					//选择一个点作为中心生成局部搜索区域
					Real rr = rnd->uniform.next();
					std::vector<Real> center;
					if (rr < 0.5) {
						center = sol1;
					}
					else {
						center = sol2;
					}
					std::vector<std::pair<Real, Real>> box;
					/*for (size_t k = 0; k < num_var; ++k) {
						std::pair<Real, Real> temp;
						temp.first = sol1[k] < sol2[k] ? sol1[k] : sol2[k];
						temp.second = sol1[k] > sol2[k] ? sol1[k] : sol2[k];
						box.emplace_back(temp);
					}*/
					for (size_t k = 0; k < num_var; ++k) {
						Real span = search_bound[k].second - search_bound[k].first;
						std::pair<Real, Real> temp;
						Real lower = center[k] - span / 40;
						temp.first = lower>search_bound[k].first?lower: search_bound[k].first;
						Real upper = center[k] + span / 40;
						temp.second = upper < search_bound[k].second ? upper : search_bound[k].second;
						box.emplace_back(temp);
					}
					m_reset_box = box;
					//auto box = getMO_HLC().subspaceTree().getBox(space);
					for (size_t j = 0; j < second_pop.size(); ++j) {
						std::vector<Real> sol;
						for (size_t j = 0; j < box.size(); ++j) {
							sol.push_back(box[j].first + rnd->uniform.next() * (box[j].second - box[j].first));
						}
						Solution<> ind4(getPop()[i][second_pop[j]]);
						ind4.variable().vect() = sol;
						ind4.evaluate(pro, alg);
						std::vector<std::shared_ptr<Solution<>>> temp_pair;
						temp_pair.emplace_back(std::make_shared<Solution<>>(getPop()[i][second_pop[j]]));
						temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
						getInteractiveSols().emplace_back(temp_pair);

						getPop()[i].getOffspring()[second_pop[j]].variable() = ind4.variable();
						getPop()[i].getOffspring()[second_pop[j]].objective() = ind4.objective();
					}
				}
				else {
					//自由种群演化
					PopDE<> temp_pop;
					std::vector<size_t> candidate, ridx;
					for (size_t j = 0; j < second_pop.size(); ++j) {
						temp_pop.append(getPop()[i][second_pop[j]]);
						candidate.push_back(j);
					}

					for (size_t j = 0; j < temp_pop.size(); ++j) {
						temp_pop.selectInCandidates(3, candidate, ridx, rnd);
						temp_pop[j].mutate(temp_pop.scalingFactor()*rnd->uniform.next(), &temp_pop[ridx[0]], &temp_pop[ridx[1]], &temp_pop[ridx[2]], pro);
						temp_pop.recombine(j, rnd, pro);

						//越界处理
						SPMOEA::repairSol(temp_pop[j].trial().variable().vect(), m_reset_box, rnd);

						temp_pop[j].trial().evaluate(pro, this);

						Solution<> ind4(temp_pop[j]);
						ind4.variable() = temp_pop[j].trial().variable();
						ind4.objective() = temp_pop[j].trial().objective();
						std::vector<std::shared_ptr<Solution<>>> temp_pair;
						temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[j]));
						temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
						getInteractiveSols().emplace_back(temp_pair);

						getPop()[i].getOffspring()[second_pop[j]].variable() = ind4.variable();
						getPop()[i].getOffspring()[second_pop[j]].objective() = ind4.objective();
					}
				}
			}
			else if (type[i] == 8) {
				size_t num_var = pro->numberVariables();
				std::vector<size_t> space_inx;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					Real v = getMO_HLC().subspaceTree().getBoxVolume(front_spaces[j]);
					if (v / getVarSpaceVolume() > std::pow(0.1, num_var)) {
						space_inx.push_back(front_spaces[j]);
					}
				}
				if (space_inx.size() > 0) {
					for (size_t j = 0; j < exploit_num; ++j) {
						size_t space_idx = (size_t)std::floor(space_inx.size() * rnd->uniform.next());
						// method 1
						/*auto link_neis = getMO_HLC().getSubspaceInfo(front_spaces[space_inx]).m_link_subspaces;
						link_neis.push_back(front_spaces[space_inx]);*/
						//method 2
						//auto space_prob=spaceLinkProbability(front_spaces[space_inx]);
						/*auto space_prob = all_front_nei_prob[space_inx];
						size_t num_neis = 1;*/

						//前沿邻域子空间交互
						auto off = sampleLinearInFrontSpace(space_inx[space_idx], j, 1, pro, alg, rnd);
						for (auto ii : off) {
							all_off.emplace_back(ii);
						}
						//统计子空间交互的成功率
						//子代在子空间的数量
						size_t sol_in_space_front = 0;
						int off_flag;//0:外面，1：前沿，-1：后排
						auto& front_sols = getMO_HLC().getSubspaceInfo(space_inx[space_idx]).m_front_sol_in_subspace;
						auto& sol = off[0];
						auto obj1 = getInteractiveSols().back().back()->objective();
						auto space_index = getMO_HLC().subspaceTree().getRegionIdx(sol);
						if (space_inx[space_idx] == space_index) {
							//判断解是不是子空间前沿
							bool flag = false;
							for (size_t p = 0; p < front_sols.size(); ++p) {
								auto& obj2 = front_sols[p]->objective();
								auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
								if (ship == Dominance::kDominated) {
									break;
								}
								else if (p == front_sols.size() - 1) {
									sol_in_space_front++;
									flag = true;
								}
							}
							if (flag) {
								off_flag = 1;
							}
							else {
								off_flag = -1;
							}
						}
						else {
							off_flag = 0;
						}
						//
						auto search = m_subspace_interactive_result.find(space_inx[space_idx]);
						if (search != m_subspace_interactive_result.end()) {
							//已有子空间
							m_subspace_interactive_result[space_inx[space_idx]].first += 1;
							m_subspace_interactive_result[space_inx[space_idx]].second.push(off_flag);
							//限制最近一定次数的子代生成结果
							while (m_subspace_interactive_result[space_inx[space_idx]].second.size() > m_queue_size) {
								m_subspace_interactive_result[space_inx[space_idx]].second.pop();
							}
						}
						else {
							//没有子空间
							std::pair<size_t, std::queue<int>> temp;
							temp.first = 1;
							temp.second.push(off_flag);
							m_subspace_interactive_result.insert(std::make_pair<>(space_inx[space_idx], temp));
						}
					}
				}
				else {
					all_off = sampleDiverse(getPop()[i], exploit_num, pro, alg, rnd);
				}
			}
			else if (type[i] == 9) {//流形学习，umap降维

			}
			if (type[i] != 6 || type[i] != 7) {
				auto& interactive_sols = getInteractiveSols();
				for (size_t k = 0; k < interactive_sols.size(); ++k) {
					auto& ind = interactive_sols[k].back();
					getPop()[i].getOffspring()[k].variable() = ind->variable();
					getPop()[i].getOffspring()[k].objective() = ind->objective();
					getPop()[i].getOffspring()[k].setCounter(0);
				}
			}
			getOperatorRatio().emplace_back(operator_ratio);
		}
		updateEE(0, num_exploit);
	}

	//基于子空间个体采样，向个体所在子空间前沿解推进，提升个体
	std::vector<std::vector<Real>> SPMOEA4_4::sampleByPush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(sol);
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable().vect() = off[0];
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//根据所在子空间位置概率选择推进方向
				auto& sol1 = sol.variable().vect();
				//将历史解分三类：支配个体，被个体支配，非支配
				std::vector<size_t> dominant_inx;
				std::vector<size_t> dominated_inx;
				std::vector<size_t> nondominated_inx;
				auto& obj1 = sol.objective();
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto& obj2 = his_sols[j]->objective();
					auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
					if (ship == Dominance::kDominant) {
						dominant_inx.push_back(j);
					}
					else if (ship == Dominance::kDominated) {
						dominated_inx.push_back(j);
					}
					else if (ship == Dominance::kNonDominated) {
						nondominated_inx.push_back(j);
					}
				}
				//根据概率选择前推或扩展
				temp_pop.append(sol);
				candidate_inx = { 0,1,2 };
				Real rand = rnd->uniform.next();
				if (rand < push_prob) {
					//前推
					if (dominated_inx.size() > 0 && dominant_inx.size() > 0) {
						//前后各加入一个
						size_t inx1 = (size_t)std::floor(dominant_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[dominant_inx[inx1]]);
						size_t inx2 = (size_t)std::floor(dominated_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[dominated_inx[inx2]]);
					}
					else if (dominant_inx.size() > 0) {
						size_t inx1 = (size_t)std::floor(dominant_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[dominant_inx[inx1]]);
						temp_pop.append(sol);
					}
					else if (dominated_inx.size() > 0) {
						temp_pop.append(sol);
						size_t inx2 = (size_t)std::floor(dominated_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[dominated_inx[inx2]]);
					}
					else {
						//加入历史解
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[inx1]);
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(*his_sols[inx2]);
					}
				}
				else {
					//扩展
					if (nondominated_inx.size() < 3) {
						//加入历史解
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[inx1]);
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(*his_sols[inx2]);
					}
					else {
						//加入非支配解
						size_t inx1 = (size_t)std::floor(nondominated_inx.size() * rnd->uniform.next());
						temp_pop.append(*his_sols[nondominated_inx[inx1]]);
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(nondominated_inx.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(*his_sols[nondominated_inx[inx2]]);
					}
				}

				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
			m_total_push_times++;
			//父代与子代的支配关系
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].trial().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_push_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_push_results[0]++;
			}
			else {
				m_push_results[2]++;
			}
		}
		return out_off;
	}

	size_t SPMOEA4_4::selectSubspace(Problem* pro, Random* rnd) {
		//综合子空间覆盖率和rank值排序，选出前排子空间
		/*std::vector<Real> volumn;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			Real v = getMO_HLC().subspaceTree().getBoxVolume(j);
			volumn.push_back(v);
		}
		Real min_v = *std::min_element(volumn.begin(), volumn.end());*/
		std::vector<size_t> split_num;
		size_t var_num = pro->numberVariables();
		size_t min_num = 1000;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//size_t num = (size_t)std::floor(min_num*volumn[j]/min_v);

			split_num.push_back(min_num);
		}
		std::vector<Real> space_coverage;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//kd-tree细分
			auto box = getMO_HLC().subspaceTree().getBox(j);
			std::vector<Real> ratios(split_num[j], 1. / split_num[j]);
			KDTree temp_tree(ratios, box);
			temp_tree.buildIndex();
			auto& his_sols = getMO_HLC().getSubspaceInfo(j).m_history_inds;
			//计算具有个体的子空间数
			std::vector<size_t> have_inds;
			for (size_t k = 0; k < his_sols.size(); ++k) {
				auto& sol = his_sols[k]->variable().vect();
				auto space_idx = temp_tree.getRegionIdx(sol);
				if (have_inds.empty()) {
					have_inds.push_back(space_idx);
				}
				else {
					if (std::find(have_inds.begin(), have_inds.end(), space_idx) == have_inds.end()) {
						have_inds.push_back(space_idx);
					}
				}
			}
			space_coverage.push_back((Real)have_inds.size() / split_num[j]);
		}
		//子空间排序值
		std::vector<Real> space_rank;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			space_rank.push_back(1 + getMO_HLC().getSubspaceInfo(j).m_best_rank);
		}
		//子空间排序
		std::vector<std::vector<Real>*> objs;
		std::vector<std::vector<Real>> all_metrics;
		for (size_t i = 0; i < space_rank.size(); ++i) {
			std::vector<Real> temp = { space_coverage[i],space_rank[i] };
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//第一层子空间
		std::vector<Real> coverage;
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(i);
				coverage.push_back(all_metrics[i][0]);
			}
		}
		//从第一层中选择覆盖率低的子空间
		auto min_cover = std::distance(coverage.begin(), std::min_element(coverage.begin(), coverage.end()));
		//size_t space = (size_t)std::floor(first_layer.size() * rnd->uniform.next());
		return first_layer[min_cover];
	}

	//根据目标空间分布，找出历史前沿中的稀疏点对
	std::map<Real,std::pair<size_t, size_t>> SPMOEA4_4::selectSubspaceFromObj(Problem* pro, Random* rnd) {
		//将目标空间所有非支配解分两类
		size_t num_obj = pro->numberObjectives();
		std::vector<std::vector<Real>> all_data;
		for (size_t i = 0; i < m_pop_index[0].size(); ++i) {
			all_data.emplace_back(getPop()[0][m_pop_index[0][i]].objective());
		}
		/*auto his_front_sols = getHisFrontSols();
		for (size_t i = 0; i < his_front_sols.size(); ++i) {
			all_data.emplace_back(his_front_sols[i]->objective());
		}*/
		auto normal_data = all_data;
		dataNormalize(normal_data);
		std::vector<Real> min_dist;
		std::vector<std::vector<Real>> all_dist;
		Real max_dist = 0;
		for (size_t i = 0; i < normal_data.size(); ++i) {
			std::vector<Real> temp_dist;
			for (size_t j = 0; j < normal_data.size(); ++j) {
				if (j != i) {
					auto dist = euclideanDistance(normal_data[i].begin(), normal_data[i].end(), normal_data[j].begin());
					temp_dist.push_back(dist);
				}
				else {
					temp_dist.push_back(0.);
				}
			}
			all_dist.emplace_back(temp_dist);
			min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
			auto temp = *std::max_element(temp_dist.begin(), temp_dist.end());
			if (temp > max_dist) {
				max_dist = temp;
			}
		}
		Real sum_dist = 0.;
		for (size_t i = 0; i < min_dist.size(); ++i) {
			sum_dist += min_dist[i];
		}
		int method = 1;
		std::vector<std::vector<size_t>> ind2cluster;//类与个体索引
		if (method == 1) {
			//dbscan聚类
			Real mean_dist = max_dist / normal_data.size();
			//将平均距离作为密度聚类的距离值
			size_t minPts = 5;//每个类中最小的个体数
			Real epsilon = 5*mean_dist;

			DBSCAN dscluster(minPts, epsilon, normal_data);
			dscluster.run();
			std::vector<int> cluster_id;//每个解所属的类,真实类和未分类等
			for (size_t i = 0; i < dscluster.m_points.size(); ++i) {
				cluster_id.push_back(dscluster.m_points[i]->clusterID);
			}
			std::vector<int> cluster_num;//所有的类
			cluster_num.push_back(cluster_id[0]);
			for (size_t i = 1; i < cluster_id.size(); ++i) {
				if (std::find(cluster_num.begin(), cluster_num.end(), cluster_id[i]) == cluster_num.end()) {
					cluster_num.push_back(cluster_id[i]);
				}
			}
			//将相同的类的点放在一起
			for (size_t i = 0; i < cluster_num.size(); ++i) {
				if (cluster_num[i] > 0) {
					std::vector<size_t> temp_cluster;
					for (size_t j = 0; j < cluster_id.size(); ++j) {
						if (cluster_id[j] == cluster_num[i]) {
							temp_cluster.push_back(j);
						}
					}
					ind2cluster.emplace_back(temp_cluster);
				}
			}
		}
		else if (method == 2) {
			//k-means聚类
			size_t K = 2;
			size_t max_iterations = 100;
			KMeans kmeans(K, normal_data.size(), num_obj, max_iterations);
			kmeans.run(normal_data, rnd);
			for (size_t i = 0; i < K; i++) {
				size_t total_points_cluster = kmeans.getClusters()[i].getTotalPoints();
				std::vector<size_t> ind_inx;
				for (size_t j = 0; j < total_points_cluster; j++) {
					ind_inx.push_back(kmeans.getClusters()[i].getPoint(j).getID());
				}
				ind2cluster.emplace_back(ind_inx);
			}
		}

		//求最近类间距离
		std::map<Real, std::pair<size_t, size_t>> output_pairs;
		if (ind2cluster.size() < 2) {
			//选稀疏点，比较每个点最近的n个点的平均距离
			size_t neigh_num = 5;
			std::vector<Real> mean_dist;
			for (size_t j = 0; j < all_dist.size(); ++j) {
				auto temp = all_dist[j];
				Real sum_d = 0;
				for (size_t k = 0; k < neigh_num; ++k) {
					Real dd = *std::min_element(temp.begin(), temp.end());
					auto inx = std::distance(temp.begin(), std::min_element(temp.begin(), temp.end()));
					sum_d += dd;
					temp[inx] = std::numeric_limits<Real>::max();
				}
				mean_dist.push_back(sum_d/neigh_num);
			}
			auto inx = std::distance(mean_dist.begin(), std::max_element(mean_dist.begin(), mean_dist.end()));
			std::pair<size_t, size_t> temp_pair;
			temp_pair.first = inx;
			temp_pair.second = inx;
			output_pairs.insert(std::make_pair<>(mean_dist[inx], temp_pair));

			////随机选两个点
			//size_t inx1 = (size_t)std::floor(ind2cluster[0].size()*rnd->uniform.next());
			//size_t inx2 = 0;
			//do {
			//	inx2 = (size_t)std::floor(ind2cluster[0].size() * rnd->uniform.next());
			//} while (inx1 == inx2);
			//Real dd = all_dist[ind2cluster[0][inx1]][ind2cluster[0][inx2]];
			//std::pair<size_t, size_t> temp_pair;
			//temp_pair.first = ind2cluster[0][inx1];
			//temp_pair.second = ind2cluster[0][inx2];
			//output_pairs.insert(std::make_pair<>(dd,temp_pair));
		}
		else {
			for (size_t i = 0; i < ind2cluster.size(); ++i) {
				auto ind_inx1 = ind2cluster[i];
				std::map <Real, std::pair<size_t, size_t>> cluster_dists;
				for (size_t j = 0; j < ind2cluster.size(); ++j) {
					if (j != i) {
						auto ind_inx2 = ind2cluster[j];
						std::vector<std::vector<Real>> temp_d;
						std::vector<Real> min_d;
						for (size_t k = 0; k < ind_inx1.size(); ++k) {
							std::vector<Real> temp;
							for (size_t p = 0; p < ind_inx2.size(); ++p) {
								temp.push_back(all_dist[ind_inx1[k]][ind_inx2[p]]);
							}
							temp_d.emplace_back(temp);
							auto min_v = *std::min_element(temp.begin(), temp.end());
							min_d.push_back(min_v);
						}
						auto mm = *std::min_element(min_d.begin(), min_d.end());
						size_t inx1 = std::distance(min_d.begin(), std::min_element(min_d.begin(), min_d.end()));
						size_t inx2 = std::distance(temp_d[inx1].begin(), std::min_element(temp_d[inx1].begin(), temp_d[inx1].end()));
						std::pair<size_t, size_t> temp_pair;
						temp_pair.first = ind_inx1[inx1];
						temp_pair.second = ind_inx2[inx2];
						cluster_dists.insert(std::make_pair<>(mm, temp_pair));
					}
				}
				//求最小距离
				auto min_pair = *std::min_element(cluster_dists.begin(), cluster_dists.end(), [](std::pair<Real, std::pair<size_t, size_t>> left, std::pair<Real, std::pair<size_t, size_t>> right) {return left.first < right.first; });
				output_pairs.insert(min_pair);
			}
		}
		
		return output_pairs;
	}

	std::vector<std::vector<Real>> SPMOEA4_4::sampleDiverse(Population<Solution<>>& pop, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_put;
		PopDE<> temp_pop;
		std::vector<size_t> candidate, ridx;
		for (size_t i = 0; i < pop.size(); ++i) {
			temp_pop.append(pop[i]);
			candidate.push_back(i);
		}

		for (size_t i = 0; i < temp_pop.size() / 2; ++i) {
			temp_pop.selectInCandidates(3, candidate, ridx, rnd);
			temp_pop[i].mutate(temp_pop.scalingFactor(), &temp_pop[ridx[0]], &temp_pop[ridx[1]], &temp_pop[ridx[2]], pro);
			temp_pop.recombine(i, rnd, pro);
			temp_pop[i].trial().evaluate(pro, this);

			Solution<> ind4(temp_pop[i]);
			ind4.variable() = temp_pop[i].trial().variable();
			ind4.objective() = temp_pop[i].trial().objective();
			out_put.emplace_back(ind4.variable().vect());
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);
		}

		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = pop.size() / 2; i < pop.size(); ++i) {
			//
			std::vector<Real> sol;
			for (size_t j = 0; j < bound.size(); ++j) {
				sol.push_back(bound[j].first + rnd->uniform.next() * (bound[j].second - bound[j].first));
			}
			Solution<> ind4(pop[i]);
			ind4.variable().vect() = sol;
			ind4.evaluate(pro, alg);
			out_put.emplace_back(sol);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(pop[i]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);
		}

		return out_put;
	}

	std::vector<std::vector<Real>> SPMOEA4_4::sampleRandom(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_put;
		//auto bound = CAST_CONOP(pro)->boundary();

		for (size_t i = 0; i < pop.size(); ++i) {
			//
			std::vector<Real> sol;
			for (size_t j = 0; j < bound.size(); ++j) {
				sol.push_back(bound[j].first + rnd->uniform.next() * (bound[j].second - bound[j].first));
			}
			Solution<> ind4(pop[i]);
			ind4.variable().vect() = sol;
			ind4.evaluate(pro, alg);
			out_put.emplace_back(sol);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(pop[i]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);
		}
		return out_put;
	}

	bool SPMOEA4_4::spaceSubdivision(Problem* pro, Random* rnd) {
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
			else if (getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity < m_divide_granularity) {
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

	void SPMOEA4_4::updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd) {
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

	void SPMOEA4_4::updateNeighSpace() {

	}

	void SPMOEA4_4::clusterSubspace() {
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

	std::vector<std::vector<size_t>> SPMOEA4_4::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA4_4::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA4_4::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	//统计子空间前沿解的交互成功率
	void SPMOEA4_4::adaptiveSplit(Problem* pro) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (auto space : m_subspace_pop_index) {
				SPMOEA_pop off_pop(0, pro);
				for (size_t k = 0; k < space.second.size(); ++k) {
					off_pop.append(getPop()[i].getOffspring()[space.second[k]]);
				}
				//统计子空间交互的成功率
				auto space_idx = space.first;
				//子代在子空间的数量
				size_t in_space = 0;
				size_t sol_in_space_front = 0;
				std::queue<int> off_flag;//0:外面，1：前沿，-1：后排
				auto& front_sols = getMO_HLC().getSubspaceInfo(space_idx).m_subspace_front_sol;
				for (size_t k = 0; k < off_pop.size(); ++k) {
					auto& sol = off_pop[k].variable().vect();
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					in_space++;
					if (space_inx == space_idx) {
						//判断解是不是子空间前沿
						bool flag = false;
						for (size_t p = 0; p < front_sols.size(); ++p) {
							auto& sol2 = front_sols[p]->variable().vect();
							if (ifSame(sol, sol2)) {
								sol_in_space_front++;
								flag = true;
								break;
							}
						}
						if (flag) {
							off_flag.push(1);
						}
						else {
							off_flag.push(-1);
						}
					}
					else {
						off_flag.push(0);
					}
				}
				//
				if (auto search = m_subspace_succ_rate.find(space_idx); search != m_subspace_succ_rate.end()) {
					//已有子空间
					m_subspace_succ_rate[space_idx].first += in_space;
					while (!off_flag.empty()) {
						m_subspace_succ_rate[space_idx].second.push(off_flag.front());
						off_flag.pop();
					}
					//限制最近的100次子代生成结果
					while (m_subspace_succ_rate[space_idx].second.size() > m_queue_size) {
						m_subspace_succ_rate[space_idx].second.pop();
					}
				}
				else {
					//没有子空间
					std::pair<size_t, std::queue<int>> temp;
					temp.first = in_space;
					temp.second = off_flag;
					m_subspace_succ_rate.insert(std::make_pair<>(space_idx, temp));
				}
			}
		}
	}

	//统计子空间前沿解的交互成功率
	void SPMOEA4_4::adaptiveSplitFrontSpace(Problem* pro, Random* rnd) {
		std::map<size_t, Real> fail_rate;
		size_t num_var = pro->numberVariables();
		for (auto& space_info : m_subspace_interactive_result) {
			if (getMO_HLC().getSubspaceInfo(space_info.first).m_best_rank == 0) {
				auto que = space_info.second.second;
				size_t fail_times = 0;
				size_t out_times = 0;
				while (!que.empty()) {
					int flag = que.front();
					if (flag == -1) {
						fail_times++;
					}
					if (flag == 0) {
						out_times++;
					}
					que.pop();
				}
				Real fail = (Real)fail_times / ((Real)space_info.second.second.size() + 0.001);
				Real out = (Real)out_times / (Real)space_info.second.second.size();
				fail_rate.insert(std::make_pair<>(space_info.first, fail));
				auto volumn = getMO_HLC().subspaceTree().getBoxVolume(space_info.first);
				if (space_info.second.second.size() >= m_queue_size && fail > 0.4 && volumn / getVarSpaceVolume() > std::pow(0.1, num_var)) {
					//先得到子空间历史解用于更新细分的子空间的信息
					int dim = findSplitDim(space_info.first, pro);
					auto& space_bound = getMO_HLC().subspaceTree().getBox(space_info.first);
					Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
					auto space_front_sols = getMO_HLC().getSubspaceInfo(space_info.first).m_subspace_front_sol;
					splitSpace(space_info.first, 2, dim, pos, false, pro, rnd);
					//根据分割后子空间前沿解的占比分配
					if (space_front_sols.size() == 0) {
						size_t a = 1;
					}
					size_t keep_num = 0;
					for (size_t j = 0; j < space_front_sols.size(); ++j) {
						if (space_front_sols[j]->variable()[dim] < pos) {
							keep_num++;
						}
					}
					keep_num = (size_t)std::floor(m_queue_size * ((Real)keep_num / (space_front_sols.size() + 0.001)));
					size_t sample_num = m_subspace_interactive_result[space_info.first].first;
					m_subspace_interactive_result[space_info.first].first = (size_t)std::floor(sample_num * keep_num / m_queue_size);
					std::queue<int> que2;
					while (m_subspace_interactive_result[space_info.first].second.size() > keep_num) {
						que2.push(m_subspace_interactive_result[space_info.first].second.front());
						m_subspace_interactive_result[space_info.first].second.pop();
					}
					auto num_spaces = getMO_HLC().numSubspace();
					std::pair<size_t, std::queue<int>> temp;
					temp.first = sample_num - m_subspace_interactive_result[space_info.first].first;
					temp.second = que2;
					m_subspace_interactive_result.insert(std::make_pair<>(num_spaces - 1, temp));
				}
			}
		}
		size_t a = 1;
	}

	void SPMOEA4_4::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA4_4::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA4_4::findSplitDim(int inx, Problem* pro) {
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
		auto bound = CAST_CONOP(pro)->boundary();
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			//dim_span.push_back(space_bound[j].second - space_bound[j].first);
			//归一化维度切分
			dim_span.push_back((space_bound[j].second - space_bound[j].first) / (bound[j].second - bound[j].first));
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

	void SPMOEA4_4::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA4_4::recordMetrics(Problem* pro, Algorithm* alg) {
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

	void SPMOEA4_4::initiObjSpace(Problem* pro) {

	}
}