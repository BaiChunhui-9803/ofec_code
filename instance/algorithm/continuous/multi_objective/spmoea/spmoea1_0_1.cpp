#include "spmoea1_0_1.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_0_1::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_0_1::run_() {
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

	void SPMOEA1_0_1::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_0_1::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(getPop().size() + 1);//最后一个为历史解
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j].phenotype());
				}
				/*for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j){
					  m_solution[1].push_back(&getPop()[i].getOffspring()[j].phenotype());
				}*/
			}
			auto& his_sols = getHisSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(&his_sols[i]->phenotype());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA1_0_1::initiVarSpace(Problem* pro) {
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

	void SPMOEA1_0_1::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//先在子空间预采样，根据比较结果生成多种群
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
		size_t num = getMO_HLC().numSubspace();
		SPMOEA_pop temp_pop(num, pro);
		std::vector<std::vector<Real>> sols;
		for (size_t i = 0; i < num; ++i) {
			auto bound = getMO_HLC().subspaceTree().getBox(i);
			for (size_t k = 0; k < bound.size(); ++k) {
				temp_pop[i].variable().vect()[k] = bound[k].first + 1. / 2 * (bound[k].second - bound[k].first);
			}
			temp_pop[i].evaluate(pro, alg);
			sols.emplace_back(temp_pop[i].variable().vect());
		}
		NDSort(temp_pop);
		//初始化子代
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
		}
		/*temp_pop.setRate(getCr(), getMr());
		temp_pop.setEta(getCeta(), getMeta());*/
		getPop().append(temp_pop);

		//更新子空间信息
		updateVarSpaceInfo(temp_pop, pro, rnd);
		SPMOEA::updateHistoryInfo(temp_pop, pro);
		SPMOEA::updateArchive(archiveNum(), pro);
		SPMOEA::updateNewPop(temp_pop);
		SPMOEA::updateObjRange(temp_pop, pro);
		updateObjSpace();
		//根据子空间排序值聚类
		clusterSubspace();

		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				auto v = getMO_HLC().subspaceTree().getBoxVolume(i);
				front_ratio += (v / total_volume);
			}
		}
		setFrontSpaceRatio(front_ratio);

		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();

		updateEE(num, 0);
		//SPMOEA::recordMetrics(pro, alg);
	}


	int SPMOEA1_0_1::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		updatePop(pro, alg, rnd);
		/********************************************************************************
				 根据各子种群的比较，设置各子种群子代的探索开发比例、子代的个体数
		********************************************************************************/
		//单纯使用子种群信息属于不完全信息比较，可否使用子空间全部的有效信息进行子空间的比较？
		//对所有子种群个体进行排序
		std::vector<Real> pop_exploit_ratio;
		std::vector<size_t> assign_pop_resource;
		PopResourceAssign(pop_exploit_ratio, assign_pop_resource, pro);
		/********************************************************************************
										   类内种群演化
		********************************************************************************/
		//1.子代中探索与开发的比例；2.子代中探索与开发的方式
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			bool m_search_balance = true;//选择子种群是否E&E平衡
			if (!m_search_balance) {
				for (auto& i : pop_exploit_ratio) {
					i = 1.;
				}
			}
			bool m_assign_resource = true;//子种群是否进行资源分配
			if (!m_search_balance) {
				for (auto& i : assign_pop_resource) {
					i = getPopsize();
				}
			}
			generateOffspring(pro, rnd, assign_pop_resource, pop_exploit_ratio);
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
		/************************************************************************
			  检测是否细分子空间:具有局部结构的子空间；位于谷地较好的子空间
		************************************************************************/
		spaceSubdivision(pro, rnd);

		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		NDSort(offspring_pop);
		SPMOEA::updateVarSpaceInfo(offspring_pop, pro, rnd);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				auto v = getMO_HLC().subspaceTree().getBoxVolume(i);
				front_ratio += (v / total_volume);
			}
		}
		setFrontSpaceRatio(front_ratio);

		/********************************************************************************************
											  种群更新，子种群内部进行淘汰选择
		*********************************************************************************************/
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			getPop()[i].envirSelection(pro, rnd, 0);
			//getPop()[i].iteration()++;

			updatePopDistribute(i, pro);
			updatePopdist(i);
		}
		/*********************************************************************************************
			 更新聚类，统计一个类中的当前个体数，补全为一个具有最小个体数的子种群，然后子种群演化
		**********************************************************************************************/
		updateCluster();

		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		return tag;
	}

	void SPMOEA1_0_1::updatePop(Problem* pro, Algorithm* alg, Random* rnd) {
		auto cur_clusters = getMO_HLC().getClusters();
		//提取原始种群个体
		std::map<size_t, std::vector<std::shared_ptr<Solution<>>>> space_inds;//每个子空间含有的当前个体
		for (size_t i = 0; i < getPop().size(); ++i) {
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
				NDSort(temp_pop0);
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
					NDSort(temp_pop1);
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
			/*temp_pop.setRate(getCr(), getMr());
			temp_pop.setEta(getCeta(), getMeta());*/
			getPop().append(temp_pop);
		}
	}

	void SPMOEA1_0_1::PopResourceAssign(std::vector<Real>& pop_exploit_ratio, std::vector<size_t>& assign_pop_resource, Problem* pro) {
		std::vector<Solution<>*> parent_pop;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); j++) {
				parent_pop.emplace_back(&getPop()[i][j]);
			}
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < parent_pop.size(); ++i) {
			objs.emplace_back(&parent_pop[i]->objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < parent_pop.size(); ++i) {
			parent_pop[i]->setFitness(rank[i]);
		}
		//根据子种群的比较情况，决定每个子种群的探索率和开发率，以及子代的个数
		std::vector<size_t> pop_min_rank;//子种群的最好rank值比较
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t min_rank = INT16_MAX;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				if (min_rank > getPop()[i][j].fitness()) {
					min_rank = getPop()[i][j].fitness();
				}
			}
			pop_min_rank.push_back(min_rank);
		}
		std::vector<Real> pop_first_ratio;//子种群前排个体占比
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t count = 0;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				if (getPop()[i][j].fitness() == pop_min_rank[i]) {
					count++;
				}
			}
			pop_first_ratio.push_back((Real)count / getPop()[i].size());
		}
		//std::vector<Real> pop_exploit_ratio;
		for (size_t i = 0; i < getPop().size(); ++i) {
			Real temp = 1 - 0.2 * pop_min_rank[i];//系数0.5决定第几层可以参与开发
			if (temp < 0) {
				temp = 0;
			}
			temp *= 0.5;//系数0.9代表上限开发率：0.9*max[0,1-0.5*min_rank]
			pop_exploit_ratio.push_back(temp);
		}
		//子空间覆盖率或者子种群的表现决定计算资源分配，即子种群产生子代的个数
		//采用C指标计算不同子种群相对其他子种群的支配率
		Matrix dominance_ratio_matrix(getPop().size(), getPop().size());
		for (size_t i = 0; i < dominance_ratio_matrix.row(); ++i) {
			for (size_t j = 0; j < dominance_ratio_matrix.col(); ++j) {
				if (i == j) {
					dominance_ratio_matrix[i][j] = 1;
				}
				else {
					//计算第j个子种群中被第i个子种群支配的个体数占j子种群个体总数的比例
					size_t count = 0;
					for (size_t p = 0; p < getPop()[j].size(); ++p) {
						bool flag = false;
						for (size_t q = 0; q < getPop()[i].size(); ++q) {
							if (getPop()[i][q].dominate(getPop()[j][p], pro)) {
								flag = true;
								break;
							}
						}
						if (flag) {
							count++;
						}
					}
					dominance_ratio_matrix[i][j] = (Real)count / getPop()[j].size();
				}
			}
		}
		//根据dominance取值范围映射计算资源
		std::vector<Real> pop_dominance;
		for (size_t i = 0; i < dominance_ratio_matrix.row(); ++i) {
			Real temp = 0.;
			for (size_t j = 0; j < dominance_ratio_matrix[i].size(); ++j) {
				temp += dominance_ratio_matrix[i][j];
			}
			pop_dominance.push_back(temp);
		}
		Real max_pop_dominance = *std::max_element(pop_dominance.begin(), pop_dominance.end());
		Real min_pop_dominance = *std::min_element(pop_dominance.begin(), pop_dominance.end());
		min_pop_dominance = 0.;
		size_t min_pop_size = 4;
		if (min_pop_size % 2 != 0) {
			min_pop_size++;
		}
		size_t max_pop_size = getPopsize();
		//std::vector<size_t> assign_pop_resource;
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t temp_num = 0;
			if (pop_dominance.size() == 1 || max_pop_dominance == min_pop_dominance) {
				temp_num = max_pop_size;
			}
			else {
				temp_num = std::floor(min_pop_size + (max_pop_size - min_pop_size) * (pop_dominance[i] - min_pop_dominance) / (max_pop_dominance - min_pop_dominance));
			}
			assign_pop_resource.push_back(temp_num);
		}
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_0_1::generateOffspring(Problem* pro, Random* rnd, const std::vector<size_t>& pop_resource, const std::vector<Real>& pop_exploit_ratio) {
		//size_t num_explore = 0;
		//size_t num_exploit = 0;
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	size_t explore_num = std::ceill(pop_resource[i] * (1 - pop_exploit_ratio[i]));
		//	size_t exploit_num = pop_resource[i] - explore_num;
		//	if (exploit_num == 0) {//使子种群至少有一个开发和探索点
		//		exploit_num++;
		//		explore_num--;
		//	}
		//	if (explore_num == 0) {//使子种群至少有一个开发和探索点
		//		explore_num++;
		//		exploit_num--;
		//	}
		//	num_explore += explore_num;
		//	num_exploit += exploit_num;
		//	//子代个体数量处理
		//	while (getPop()[i].getOffspring().size() - getPop()[i].size() > pop_resource[i]) {
		//		getPop()[i].getOffspring().remove(getPop()[i].getOffspring().end() - 1);
		//	}
		//	while (getPop()[i].getOffspring().size() - getPop()[i].size() < pop_resource[i]) {
		//		Solution<> temp_ind(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
		//		getPop()[i].getOffspring().append(temp_ind);
		//	}
		//	auto spaces = getMO_HLC().getCluster(i);
		//	std::vector<size_t> total_spaces;
		//	for (size_t j = 0; j < spaces.size(); ++j) {
		//		for (size_t k = 0; k < spaces[j].size(); ++k) {
		//			total_spaces.push_back(spaces[j][k]);
		//		}
		//	}
		//	//在子种群前排个体间开发,开发方式选择与邻域个体之间交互
		//	auto search_bound = CAST_CONOP(pro)->boundary();
		//	bool neighbor_operator = false;
		//	size_t front_num = 0;
		//	for (size_t j = 0; j < spaces.size(); ++j) {
		//		for (size_t k = 0; k < spaces[j].size(); ++k) {
		//			auto& front_sols = getMO_HLC().getSubspaceInfo(spaces[j][k]).m_subspace_front_sol;
		//			for (size_t p = 0; p < front_sols.size(); ++p) {
		//				front_num++;
		//			}
		//		}
		//	}
		//	SPMOEA_pop front_pop(front_num, pro);
		//	size_t count_num = 0;
		//	for (size_t j = 0; j < spaces.size(); ++j) {
		//		for (size_t k = 0; k < spaces[j].size(); ++k) {
		//			auto& front_sols = getMO_HLC().getSubspaceInfo(spaces[j][k]).m_subspace_front_sol;
		//			for (size_t p = 0; p < front_sols.size(); ++p) {
		//				front_pop[count_num] = *getMO_HLC().getSubspaceInfo(spaces[j][k]).m_subspace_front_sol[p];
		//				count_num++;
		//			}
		//		}
		//	}
		//	for (size_t j = 0; j < exploit_num; ++j) {
		//		if (neighbor_operator && front_num > 4) {
		//			//选择前沿子空间与其邻域子空间内的个体交互
		//			front_pop.setRate(1, 1. / CAST_CONOP(pro)->numberVariables());
		//			front_pop.setEta(10, 10);
		//			std::vector<size_t> p(2);
		//			do {//不重复采样
		//				do {
		//					p[0] = front_pop.tournamentSelection(pro, rnd);
		//					p[1] = front_pop.tournamentSelection(pro, rnd);
		//				} while (p[1] == p[0]);
		//				front_pop.crossover(p[0], p[1], front_pop.getOffspring()[j], front_pop.getOffspring()[j + 1], pro, rnd);
		//				front_pop.mutate(front_pop.getOffspring()[j], pro, rnd);
		//				front_pop.mutate(front_pop.getOffspring()[j + 1], pro, rnd);
		//				//子代越界处理
		//				repairSol(front_pop.getOffspring()[j].variable().vect(), total_spaces, rnd);
		//				repairSol(front_pop.getOffspring()[j + 1].variable().vect(), total_spaces, rnd);
		//			} while (front_pop.getOffspring()[j].variable().vect() == front_pop[p[0]].variable().vect() || \
		//				front_pop.getOffspring()[j].variable().vect() == front_pop[p[1]].variable().vect() || \
		//				front_pop.getOffspring()[j + 1].variable().vect() == front_pop[p[0]].variable().vect() || \
		//				front_pop.getOffspring()[j + 1].variable().vect() == front_pop[p[1]].variable().vect());

		//			front_pop.getOffspring()[j].setCounter(0);
		//			front_pop.getOffspring()[j + 1].setCounter(0);

		//			//每次确定一个解，从2个子代中随机选择一个
		//			Real rnd_num = rnd->uniform.next();
		//			if (rnd_num > 0.5) {
		//				getPop()[i].getOffspring()[j].variable() = front_pop.getOffspring()[j + 1].variable();
		//			}
		//		}
		//		else {
		//			std::vector<size_t> p(2);
		//			do {//不重复采样
		//				do {
		//					p[0] = getPop()[i].tournamentSelection(pro, rnd);
		//					p[1] = getPop()[i].tournamentSelection(pro, rnd);
		//				} while (p[1] == p[0]);
		//				getPop()[i].crossover(p[0], p[1], getPop()[i].getOffspring()[j], getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				//子代越界处理
		//				repairSol(getPop()[i].getOffspring()[j].variable().vect(), total_spaces, rnd);
		//				repairSol(getPop()[i].getOffspring()[j + 1].variable().vect(), total_spaces, rnd);
		//			} while (getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[1]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[1]].variable().vect());

		//			getPop()[i].getOffspring()[j].setCounter(0);
		//			getPop()[i].getOffspring()[j + 1].setCounter(0);

		//			//每次确定一个解，从2个子代中随机选择一个
		//			Real rnd_num = rnd->uniform.next();
		//			if (rnd_num > 0.5) {
		//				getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//			}
		//		}

		//	}
		//	//在子种群所在吸引域内的子空间内进行探索
		//	//均匀探索，还是计算子空间的覆盖率
		//	bool random_sample = false;
		//	std::vector<size_t> sample_spaces;
		//	if (random_sample) {//随机采样
		//		for (size_t j = 0; j < explore_num; ++j) {
		//			size_t ind_inx = std::floor(total_spaces.size() * rnd->uniform.next());
		//			sample_spaces.push_back(total_spaces[ind_inx]);
		//		}
		//	}
		//	else {//基于覆盖率采样
		//		//确定每个子空间的划分份数，考虑子空间的大小均匀性和子空间内历史解的个数
		//		std::vector<Real> space_density;
		//		for (size_t j = 0; j < total_spaces.size(); ++j) {
		//			auto volume = getMO_HLC().subspaceTree().getBoxVolume(total_spaces[j]);
		//			space_density.push_back(getMO_HLC().getSubspaceInfo(total_spaces[j]).m_history_inds.size() / volume);
		//		}
		//		//在密度最小的子空间采样
		//		size_t select_count = 0;
		//		while (select_count < explore_num) {
		//			//找最低的覆盖率的子空间
		//			size_t min_inx = std::distance(space_density.begin(), std::min_element(space_density.begin(), space_density.end()));
		//			sample_spaces.push_back(total_spaces[min_inx]);
		//			//更新子空间的覆盖率
		//			auto volume = getMO_HLC().subspaceTree().getBoxVolume(total_spaces[min_inx]);
		//			space_density[min_inx] = space_density[min_inx] + 1. / volume;
		//			select_count++;
		//		}

		//		//size_t max_inx = std::distance(space_density.begin(),std::max_element(space_density.begin(), space_density.end()));
		//		//size_t base_num = getMO_HLC().getSubspaceInfo(total_spaces[max_inx]).m_history_inds.size();
		//		//Real unit_volume = getMO_HLC().subspaceTree().getBoxVolume(total_spaces[max_inx])/base_num;
		//		//std::vector<Real> space_split_num;
		//		//for (size_t j = 0; j < spaces.size(); ++j) {
		//		//	if (j == max_inx) {
		//		//		space_split_num.push_back(base_num);
		//		//	}
		//		//	else {
		//		//		auto num_sp = std::ceil(getMO_HLC().subspaceTree().getBoxVolume(total_spaces[j]) / unit_volume);
		//		//		space_split_num.push_back(num_sp);
		//		//	}
		//		//}
		//		//std::vector<Real> space_cover_rate;
		//		//for (size_t j = 0; j < spaces.size(); ++j) {
		//		//	space_cover_rate.push_back(spaceCoverRatio(total_spaces[j],space_split_num[j]));
		//		//}
		//		//size_t select_count = 0;
		//		//while (select_count < explore_num) {
		//		//	//找最低的覆盖率的子空间
		//		//	auto sp_inx = std::distance(space_cover_rate.begin(),std::min_element(space_cover_rate.begin(),space_cover_rate.end()));
		//		//	sample_spaces.push_back(total_spaces[sp_inx]);
		//		//	//更新子空间的覆盖率
		//		//	space_cover_rate[sp_inx] = (space_split_num[sp_inx] * space_cover_rate[sp_inx] + 1) / space_split_num[sp_inx];
		//		//	select_count++;
		//		//}
		//	}

		//	//在子空间均匀随机采样还是基于已有个体产生新解
		//	bool sample_random = false;
		//	if (sample_random) {
		//		for (size_t j = 0; j < sample_spaces.size(); ++j) {
		//			auto bound = getMO_HLC().subspaceTree().getBox(sample_spaces[j]);
		//			for (size_t k = 0; k < bound.size(); ++k) {
		//				getPop()[i].getOffspring()[exploit_num + j].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
		//			}
		//		}
		//	}
		//	else {
		//		for (size_t j = 0; j < sample_spaces.size(); ++j) {
		//			auto bound = getMO_HLC().subspaceTree().getBox(sample_spaces[j]);
		//			auto& front_sols = getMO_HLC().getSubspaceInfo(sample_spaces[j]).m_subspace_front_sol;
		//			auto& his_sols = getMO_HLC().getSubspaceInfo(sample_spaces[j]).m_history_inds;
		//			if (front_sols.size() < 2 || his_sols.size() < 4) {
		//				//随机采样一个解
		//				for (size_t k = 0; k < bound.size(); ++k) {
		//					getPop()[i].getOffspring()[exploit_num + j].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
		//				}
		//			}
		//			else {
		//				SPMOEA_pop temp_pop(2, pro);
		//				temp_pop.setRate(1, 1. / bound.size());
		//				temp_pop.setEta(10, 10);
		//				std::vector<Real> sol1, sol2;
		//				size_t front_inx = 0;
		//				size_t his_inx = 0;
		//				do {
		//					do {
		//						front_inx = std::floor(front_sols.size() * rnd->uniform.next());
		//						his_inx = std::floor(his_sols.size() * rnd->uniform.next());
		//						sol1 = front_sols[front_inx]->variable().vect();
		//						sol2 = his_sols[his_inx]->variable().vect();
		//					} while (sol1 == sol2);
		//					temp_pop[0] = *front_sols[front_inx];
		//					temp_pop[1] = *his_sols[his_inx];
		//					temp_pop.crossover(0, 1, temp_pop.getOffspring()[0], temp_pop.getOffspring()[1], pro, rnd);
		//					temp_pop.mutate(temp_pop.getOffspring()[0], pro, rnd);
		//					temp_pop.mutate(temp_pop.getOffspring()[1], pro, rnd);
		//					//子代越界处理
		//					repairSol(temp_pop.getOffspring()[0].variable().vect(), bound, rnd);
		//					repairSol(temp_pop.getOffspring()[1].variable().vect(), bound, rnd);
		//				} while (temp_pop.getOffspring()[0].variable().vect() == temp_pop[0].variable().vect() || \
		//					temp_pop.getOffspring()[0].variable().vect() == temp_pop[1].variable().vect() || \
		//					temp_pop.getOffspring()[1].variable().vect() == temp_pop[0].variable().vect() || \
		//					temp_pop.getOffspring()[1].variable().vect() == temp_pop[1].variable().vect());

		//				//每次确定一个解，从2个子代中随机选择一个
		//				Real rnd_num = rnd->uniform.next();
		//				if (rnd_num > 0.5) {
		//					getPop()[i].getOffspring()[exploit_num + j].variable() = temp_pop.getOffspring()[0].variable();
		//				}
		//				else {
		//					getPop()[i].getOffspring()[exploit_num + j].variable() = temp_pop.getOffspring()[1].variable();
		//				}
		//			}
		//		}
		//	}
		//}
		//updateEE(num_explore, num_exploit);
	}

	void SPMOEA1_0_1::spaceSubdivision(Problem* pro, Random* rnd) {
		Real total_volume = getVarSpaceVolume();
		size_t num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		//细分边界子空间
		bool split_boundary = 0;
		for (size_t i = 0; i < num_spaces; ++i) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(i);
			auto belong_clusters = getMO_HLC().getSubspaceInfo(i).idx_cluster;
			if (split_boundary) {
				if (belong_clusters.size() > 1) {
					//根据子空间大小决定是否细分
					if (space_volume / total_volume > std::pow(0.1, num_var)) {
						//先得到子空间历史解用于更新细分的子空间的信息
						splitSpace(i, 2 * belong_clusters.size(), 0, 0, true, pro, rnd);
					}
				}
			}
		}
		//细分全局前沿子空间还是子种群前沿子空间
		bool split_front = 0;
		auto pre_clusters = getMO_HLC().getClusters();
		if (split_front) {
			for (size_t i = 0; i < num_spaces; ++i) {
				auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(i);
				if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
					if (space_volume / total_volume > std::pow(0.05, num_var)) {
						//先得到子空间历史解用于更新细分的子空间的信息
						int dim = findSplitDim(i, pro);
						auto& space_bound = getMO_HLC().subspaceTree().getBox(i);
						Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
						splitSpace(i, 2, dim, pos, false, pro, rnd);
					}
				}
			}
		}
		else {
			for (size_t i = 0; i < pre_clusters.size(); ++i) {
				auto front_spaces = pre_clusters[i][0];
				//int temp_rank = getMO_HLC().getSubspaceInfo(pre_clusters[i][0][0]).m_best_rank;
				for (size_t j = 0; j < front_spaces.size(); ++j) {
					auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(front_spaces[j]);
					Real min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
					min_num = 10;
					if (space_volume / total_volume > std::pow(1. / min_num, num_var)) {
						//先得到子空间历史解用于更新细分的子空间的信息
						int dim = findSplitDim(front_spaces[j], pro);
						auto& space_bound = getMO_HLC().subspaceTree().getBox(front_spaces[j]);
						Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
						splitSpace(front_spaces[j], 2, dim, pos, false, pro, rnd);
					}
				}
			}
		}
		//在局部结构中寻找新的局部结构
		bool find_new_region = 0;
		//根据种群比较的结果，检测位于前沿的种群所在的子空间是否具有可分性
		if (find_new_region) {
			bool in_optima_space = true;//在前沿子空间寻找局部结构
			if (in_optima_space) {
				for (size_t i = 0; i < pre_clusters.size(); ++i) {
					if (getMO_HLC().getSubspaceInfo(pre_clusters[i][0][0]).m_best_rank == 0) {
						//先得到子空间历史解用于更新细分的子空间的信息
						for (size_t j = 0; j < pre_clusters[i][0].size(); ++j) {
							auto split_flag = subspaceSeperable(pre_clusters[i][0][j], pro);//聚类子空间是否可分
							if (split_flag > 0) {
								splitSpace(pre_clusters[i][0][j], 8, 0, 0, true, pro, rnd);//细分2次幂
							}
						}
					}
					else {
						break;
					}
				}
			}
			else {
				for (size_t i = 0; i < pre_clusters.size(); ++i) {
					for (size_t j = 0; j < pre_clusters[i][0].size(); ++j) {
						auto split_flag = subspaceSeperable(pre_clusters[i][0][j], pro);//聚类子空间是否可分
						if (split_flag > 0) {
							splitSpace(pre_clusters[i][0][j], 8, 0, 0, true, pro, rnd);//细分2次幂
						}
					}
				}
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

	void SPMOEA1_0_1::updateCluster() {
		clusterSubspace();
		//逐个比较每个子空间及其邻域的好坏，发现局部结构
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		//先找出排序值比周围子空间好的子空间
		std::vector<size_t> local_good_space;
		bool strict_nei = false;
		for (size_t i = 0; i < space_num; ++i) {
			size_t cur_rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
			auto neigh_inx = getMO_HLC().getSubspaceInfo(i).m_sub_neighbors;
			std::list<int> neighs;
			if (strict_nei) {
				//满足一定重合率的子空间邻域
				for (auto k : neigh_inx) {
					auto center_box = getMO_HLC().subspaceTree().getBox(i);
					auto compara_box = getMO_HLC().subspaceTree().getBox(k);
					size_t dim = INT16_MAX;
					for (size_t p = 0; p < center_box.size(); ++p) {
						auto bound1 = center_box[p];
						auto bound2 = compara_box[p];
						//判断该维是否重叠
						if (bound1.first < bound2.second && bound1.second > bound2.first) {
							dim = p;
							break;
						}
					}
					if (dim < INT16_MAX) {
						neighs.push_back(k);
					}
				}
			}
			else {
				for (auto k : neigh_inx) {
					neighs.push_back(k);
				}
			}
			bool is_local = true;
			for (auto j : neighs) {
				if (getMO_HLC().getSubspaceInfo(j).m_best_rank < cur_rank) {
					is_local = false;
					break;
				}
			}
			if (is_local) {
				local_good_space.push_back(i);
			}
		}
		//对局部前沿子空间聚类
		std::vector<size_t> flag_idx(space_num, 0);
		std::vector<std::vector<size_t>> cluster_idx;
		cluster_idx = clusterFrontSpace(local_good_space);
		for (size_t i = 0; i < cluster_idx.size(); ++i) {
			for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
				flag_idx[cluster_idx[i][j]] = 1;
			}
		}
		//保存当前每个类的最外层的最大的rank值
		std::vector<int> cluster_cur_ranks;
		for (size_t i = 0; i < cluster_idx.size(); ++i) {
			auto temp_rank = getMO_HLC().getSubspaceInfo(cluster_idx[i][0]).m_best_rank;
			cluster_cur_ranks.push_back(temp_rank);
		}
		std::vector<std::vector<std::vector<size_t>>> final_clusters;
		for (size_t i = 0; i < cluster_idx.size(); ++i) {
			std::vector<std::vector<size_t>> temp_cluster;
			temp_cluster.emplace_back(cluster_idx[i]);
			final_clusters.emplace_back(temp_cluster);
		}
		//根据每个类的最外层rank值的大小决定哪个类向外扩散一层
		while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
			//先找出具有最小最大rank值的类的索引
			auto itr = std::min_element(cluster_cur_ranks.begin(), cluster_cur_ranks.end());
			auto cluster_inx = std::distance(cluster_cur_ranks.begin(), itr);
			//当前类的所有子空间索引
			std::vector<size_t> current_cluster_inx;
			for (size_t i = 0; i < final_clusters[cluster_inx].size(); ++i) {
				for (size_t j = 0; j < final_clusters[cluster_inx][i].size(); ++j) {
					current_cluster_inx.push_back(final_clusters[cluster_inx][i][j]);
				}
			}
			//然后找出它们邻域的索引
			auto cur_cluster = final_clusters[cluster_inx].back();//当前比较层
			std::vector<size_t> neigh_inx;//当前比较层的外层邻域
			for (size_t i = 0; i < cur_cluster.size(); ++i) {
				auto temp_nei = getMO_HLC().getSubspaceInfo(cur_cluster[i]).m_sub_neighbors;
				std::list<int> neighs;
				if (strict_nei) {
					//满足一定重合率的子空间邻域
					for (auto k : temp_nei) {
						auto center_box = getMO_HLC().subspaceTree().getBox(cur_cluster[i]);
						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
						size_t dim = INT16_MAX;
						for (size_t p = 0; p < center_box.size(); ++p) {
							auto bound1 = center_box[p];
							auto bound2 = compara_box[p];
							//判断该维是否重叠
							if (bound1.first < bound2.second && bound1.second > bound2.first) {
								dim = p;
								break;
							}
						}
						if (dim < INT16_MAX) {
							neighs.push_back(k);
						}
					}
				}
				else {
					for (auto k : temp_nei) {
						neighs.push_back(k);
					}
				}
				for (auto j : neighs) {
					if (std::find(current_cluster_inx.begin(), current_cluster_inx.end(), j) == current_cluster_inx.end()) {
						if (std::find(neigh_inx.begin(), neigh_inx.end(), j) == neigh_inx.end()) {
							neigh_inx.push_back(j);
						}
					}
				}
			}
			//再判断邻域索引子空间是否被当前层的支配，若是，加入聚类，若不是，不加入
			bool permit_overlap = true; //是否允许重叠聚类子空间
			std::vector<size_t> temp_layer;
			for (size_t i = 0; i < neigh_inx.size(); ++i) {
				bool add_in = false;
				if (permit_overlap || flag_idx[neigh_inx[i]] == 0) {
					add_in = true;
				}
				else {
					add_in = false;
				}
				if (add_in) {
					size_t cur_rank = getMO_HLC().getSubspaceInfo(neigh_inx[i]).m_best_rank;
					std::vector<size_t> neigh_cur_compara;
					auto nei_inx = getMO_HLC().getSubspaceInfo(neigh_inx[i]).m_sub_neighbors;
					std::list<int> neighs;
					if (strict_nei) {
						//满足一定重合率的子空间邻域
						for (auto k : nei_inx) {
							auto center_box = getMO_HLC().subspaceTree().getBox(neigh_inx[i]);
							auto compara_box = getMO_HLC().subspaceTree().getBox(k);
							size_t dim = INT16_MAX;
							for (size_t p = 0; p < center_box.size(); ++p) {
								auto bound1 = center_box[p];
								auto bound2 = compara_box[p];
								//判断该维是否重叠
								if (bound1.first < bound2.second && bound1.second > bound2.first) {
									dim = p;
									break;
								}
							}
							if (dim < INT16_MAX) {
								neighs.push_back(k);
							}
						}
					}
					else {
						for (auto k : nei_inx) {
							neighs.push_back(k);
						}
					}
					for (size_t j = 0; j < cur_cluster.size(); ++j) {
						if (std::find(neighs.begin(), neighs.end(), cur_cluster[j]) != neighs.end()) {
							neigh_cur_compara.push_back(cur_cluster[j]);
						}
					}
					for (size_t j = 0; j < neigh_cur_compara.size(); ++j) {//邻域子空间rank值比较
						size_t compara_rank = getMO_HLC().getSubspaceInfo(neigh_cur_compara[j]).m_best_rank;
						if (cur_rank >= compara_rank) {
							temp_layer.push_back(neigh_inx[i]);
							flag_idx[neigh_inx[i]] = 1;
							break;
						}
					}
				}
			}
			if (!temp_layer.empty()) {
				final_clusters[cluster_inx].emplace_back(temp_layer);
			}
			//最后更新当前类的最外层的rank
			int cur_layer_rank = -1 * INT16_MAX;
			if (temp_layer.empty()) {
				cluster_cur_ranks[cluster_inx] = INT16_MAX;
			}
			else {
				for (size_t i = 0; i < temp_layer.size(); ++i) {
					int temp_rank = getMO_HLC().getSubspaceInfo(temp_layer[i]).m_best_rank;
					if (temp_rank > cur_layer_rank) {
						cur_layer_rank = temp_rank;
					}
				}
				cluster_cur_ranks[cluster_inx] = cur_layer_rank;
			}
		}
		//查看每个类的第一层是否被其他类包含，若是，则合并
		size_t num_eval = m_evaluations;
		m_his_effect_eval.push_back(num_eval);
		std::vector<bool> effective_flag(final_clusters.size(), true);
		for (size_t i = 0; i < final_clusters.size(); ++i) {
			auto front_layer = final_clusters[i][0][0];
			for (size_t j = 0; j < final_clusters.size(); ++j) {
				if (j != i) {
					bool flag = false;
					for (size_t k = 0; k < final_clusters[j].size(); ++k) {
						if (std::find(final_clusters[j][k].begin(), final_clusters[j][k].end(), front_layer) != final_clusters[j][k].end()) {
							effective_flag[i] = false;
							flag = true;
							break;
						}
					}
					if (flag) {
						break;
					}
				}
			}
		}
		std::vector<std::vector<std::vector<size_t>>> final_effect_clusters;
		for (size_t i = 0; i < final_clusters.size(); ++i) {
			if (effective_flag[i]) {
				final_effect_clusters.emplace_back(final_clusters[i]);
			}
		}
		for (size_t i = 0; i < final_effect_clusters.size(); ++i) {
			getMO_HLC().getClusters().push_back(final_effect_clusters[i]);
		}
		//更新子空间所属类
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
		}
		for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					getMO_HLC().getSubspaceInfo(getMO_HLC().getClusters()[i][j][k]).idx_cluster.push_back(i);
				}
			}
		}
	}

	//void SPMOEA1_0_1::clusterSubspace() {//最初版本
	//	//每次聚类只加入一层，允许类间有重叠
	//	getMO_HLC().getClusters().clear();
	//	auto space_num = getMO_HLC().numSubspace();
	//	std::vector<size_t> flag_idx(space_num, 0);
	//	//先找出排序值最靠前的子空间
	//	while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
	//		//先得到现有子空间最好的排序值的索引
	//		int rank = INT16_MAX;//存放一个类最靠前的rank值
	//		std::vector<size_t> best_idx;//当前聚类具有最好rank值的子空间索引
	//		for (size_t i = 0; i < space_num; ++i) {
	//			if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank < rank) {
	//				rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
	//			}
	//		}
	//		for (size_t i = 0; i < space_num; ++i) {
	//			if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank == rank)
	//				best_idx.push_back(i);
	//		}
	//		//对这些最靠前的子空间进行聚类
	//		std::vector<std::vector<size_t>> cluster_idx;
	//		cluster_idx = clusterFrontSpace(best_idx);

	//		for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//			for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
	//				flag_idx[cluster_idx[i][j]] = 1;
	//			}
	//		}

	//		bool gradual_cluster = true;//是否逐层扩散
	//		bool strict_cluster = true;//对于每个子空间，只有其未聚类的邻域排序值全部小于该子空间，才加入到聚类
	//		bool permit_overlap = true;//是否允许聚类的子空间重叠
	//		if (gradual_cluster) {
	//			auto current = cluster_idx;
	//			while (!current.empty()) {
	//				std::vector<std::vector<size_t>> temp_current;//所有类的当前比较的子空间
	//				for (size_t i = 0; i < current.size(); ++i) {
	//					std::vector<size_t> temp1;//每一个类的第一层邻域
	//					if (current[i].size() > 1 || current[i][0] != INT16_MAX) {
	//						for (size_t j = 0; j < current[i].size(); ++j) {
	//							size_t temp_rank = getMO_HLC().getSubspaceInfo(current[i][j]).m_best_rank;
	//							auto temp_neighbor = getMO_HLC().getSubspaceInfo(current[i][j]).m_sub_neighbors;
	//							std::vector<size_t> sub;
	//							for (auto& k : temp_neighbor) {
	//								if (permit_overlap) {//是否允许聚类的子空间重叠
	//									sub.push_back(k);
	//								}
	//								else {
	//									if (flag_idx[k] == 0) {
	//										sub.push_back(k);
	//									}
	//								}
	//							}
	//							if (strict_cluster) {
	//								std::vector<Real> temp_space;
	//								bool add = true;
	//								for (size_t p = 0; p < sub.size(); ++p) {
	//									if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
	//										if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank > temp_rank) {
	//											temp_space.push_back(sub[p]);
	//										}
	//										else {
	//											add = false;
	//											break;
	//										}
	//									}
	//								}
	//								if (add) {
	//									for (auto& k : temp_space) {
	//										cluster_idx[i].push_back(k);
	//										temp1.push_back(k);
	//										flag_idx[k] = 1;
	//									}
	//								}
	//							}
	//							else {
	//								for (size_t p = 0; p < sub.size(); ++p) {
	//									if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank > temp_rank) {
	//										if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
	//											cluster_idx[i].push_back(sub[p]);
	//											temp1.push_back(sub[p]);
	//											flag_idx[sub[p]] = 1;
	//										}
	//									}
	//								}
	//							}
	//						}
	//					}
	//					if (temp1.empty()) {
	//						temp1.push_back(INT16_MAX);
	//					}
	//					temp_current.push_back(temp1);
	//				}
	//				bool over = false;
	//				for (size_t i = 0; i < temp_current.size(); ++i) {
	//					if (temp_current[i][0] != INT16_MAX) {
	//						break;
	//					}
	//					else if (i == temp_current.size() - 1) {
	//						over = true;
	//					}
	//				}
	//				if (over) {
	//					current.clear();
	//				}
	//				else {
	//					current = temp_current;
	//					temp_current.clear();
	//				}
	//			}
	//		}
	//		else {
	//			for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//				std::vector<size_t> cur_compare_space = cluster_idx[i];
	//				while (!empty(cur_compare_space)) {
	//					auto compare_space = cur_compare_space;
	//					cur_compare_space.clear();
	//					for (size_t j = 0; j < compare_space.size(); ++j) {
	//						auto neighbors = getMO_HLC().getSubspaceInfo(compare_space[j]).m_sub_neighbors;
	//						std::vector<size_t> sub;
	//						for (auto& k : neighbors) {
	//							if (permit_overlap) {//是否允许聚类的子空间重叠
	//								sub.push_back(k);
	//							}
	//							else {
	//								if (flag_idx[k] == 0) {
	//									sub.push_back(k);
	//								}
	//							}
	//						}
	//						if (strict_cluster) {
	//							std::vector<Real> temp_space;
	//							bool add = true;
	//							for (size_t p = 0; p < sub.size(); ++p) {
	//								if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
	//									if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank > getMO_HLC().getSubspaceInfo(compare_space[j]).m_best_rank) {
	//										temp_space.push_back(sub[p]);
	//										cluster_idx[i].push_back(sub[p]);
	//										cur_compare_space.push_back(sub[p]);
	//										flag_idx[sub[p]] = 1;
	//									}
	//									else {
	//										add = false;
	//										break;
	//									}
	//								}
	//							}
	//							if (add) {
	//								for (auto& k : temp_space) {
	//									cluster_idx[i].push_back(k);
	//									cur_compare_space.push_back(k);
	//									flag_idx[k] = 1;
	//								}
	//							}
	//						}
	//						else {
	//							for (size_t p = 0; p < sub.size(); ++p) {
	//								if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank > getMO_HLC().getSubspaceInfo(compare_space[j]).m_best_rank) {
	//									if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
	//										cluster_idx[i].push_back(sub[p]);
	//										cur_compare_space.push_back(sub[p]);
	//										flag_idx[sub[p]] = 1;
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//		for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//			getMO_HLC().getClusters().push_back(cluster_idx[i]);
	//		}
	//	}
	//	//更新子空间所属类
	//	for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
	//		getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
	//	}
	//	for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
	//		for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
	//			getMO_HLC().getSubspaceInfo(getMO_HLC().getClusters()[i][j]).idx_cluster.push_back(i);
	//		}
	//	}
	//}

 //   void SPMOEA1_0_1::clusterSubspace() {
	//	//每次聚类只加入一层，允许类间有重叠
	//	getMO_HLC().getClusters().clear();
	//	auto space_num = getMO_HLC().numSubspace();
	//	//先找出排序值比周围子空间好的子空间
	//	std::vector<size_t> local_good_space;
	//	bool strict_nei = false;
	//	for (size_t i = 0; i < space_num; ++i) {
	//		size_t cur_rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
	//		auto neigh_inx = getMO_HLC().getSubspaceInfo(i).m_sub_neighbors;
	//		std::list<int> neighs;
	//		if (strict_nei) {
	//			//满足一定重合率的子空间邻域
	//			for (auto k : neigh_inx) {
	//				auto center_box = getMO_HLC().subspaceTree().getBox(i);
	//				auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//				size_t dim = INT16_MAX;
	//				for (size_t p = 0; p < center_box.size(); ++p) {
	//					auto bound1 = center_box[p];
	//					auto bound2 = compara_box[p];
	//					//判断该维是否重叠
	//					if (bound1.first < bound2.second && bound1.second > bound2.first) {
	//						dim = p;
	//						break;
	//					}
	//				}
	//				if (dim < INT16_MAX) {
	//					neighs.push_back(k);
	//				}
	//			}
	//		}
	//		else {
	//			for (auto k : neigh_inx) {
	//				neighs.push_back(k);
	//			}
	//		}
	//		bool is_local=true;
	//		for (auto j : neighs) {
	//			if (getMO_HLC().getSubspaceInfo(j).m_best_rank < cur_rank) {
	//				is_local = false;
	//				break;
	//			}
	//		}
	//		if (is_local) {
	//			local_good_space.push_back(i);
	//		}
	//	}
	//	//对局部前沿子空间聚类
	//	std::vector<size_t> flag_idx(space_num, 0);
	//	std::vector<std::vector<size_t>> cluster_idx;
	//	cluster_idx = clusterFrontSpace(local_good_space);
	//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//		for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
	//			flag_idx[cluster_idx[i][j]] = 1;
	//		}
	//	}
	//	//保存当前每个类的最外层的最大的rank值
	//	std::vector<int> cluster_cur_ranks;
	//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//		auto temp_rank= getMO_HLC().getSubspaceInfo(cluster_idx[i][0]).m_best_rank;
	//		cluster_cur_ranks.push_back(temp_rank);
	//	}
	//	std::vector<std::vector<std::vector<size_t>>> final_clusters;
	//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//		std::vector<std::vector<size_t>> temp_cluster;
	//		temp_cluster.emplace_back(cluster_idx[i]);
	//		final_clusters.emplace_back(temp_cluster);
	//	}
	//	//根据每个类的最外层rank值的大小决定哪个类向外扩散一层
	//	while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
	//		//先找出具有最小最大rank值的类的索引
	//		auto itr=std::min_element(cluster_cur_ranks.begin(),cluster_cur_ranks.end());
	//		auto cluster_inx = std::distance(cluster_cur_ranks.begin(),itr);
	//		//当前类的所有子空间索引
	//		std::vector<size_t> current_cluster_inx;
	//		for (size_t i = 0; i < final_clusters[cluster_inx].size(); ++i) {
	//			for (size_t j = 0; j < final_clusters[cluster_inx][i].size(); ++j) {
	//				current_cluster_inx.push_back(final_clusters[cluster_inx][i][j]);
	//			}
	//		}
	//		//然后找出它们邻域的索引
	//		auto cur_cluster = final_clusters[cluster_inx].back();//当前比较层
	//		std::vector<size_t> neigh_inx;//当前比较层的外层邻域
	//		for (size_t i = 0; i < cur_cluster.size(); ++i) {
	//			auto temp_nei = getMO_HLC().getSubspaceInfo(cur_cluster[i]).m_sub_neighbors;
	//			std::list<int> neighs;
	//			if (strict_nei) {
	//				//满足一定重合率的子空间邻域
	//				for (auto k : temp_nei) {
	//					auto center_box = getMO_HLC().subspaceTree().getBox(cur_cluster[i]);
	//					auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//					size_t dim = INT16_MAX;
	//					for (size_t p = 0; p < center_box.size(); ++p) {
	//						auto bound1 = center_box[p];
	//						auto bound2 = compara_box[p];
	//						//判断该维是否重叠
	//						if (bound1.first < bound2.second && bound1.second > bound2.first) {
	//							dim = p;
	//							break;
	//						}
	//					}
	//					if (dim < INT16_MAX) {
	//						neighs.push_back(k);
	//					}
	//				}
	//			}
	//			else {
	//				for (auto k : temp_nei) {
	//					neighs.push_back(k);
	//				}
	//			}
	//			for (auto j : neighs) {
	//				if (std::find(current_cluster_inx.begin(), current_cluster_inx.end(), j) == current_cluster_inx.end()) {
	//					if (std::find(neigh_inx.begin(), neigh_inx.end(), j) == neigh_inx.end()) {
	//						neigh_inx.push_back(j);
	//					}
	//				}
	//			}
	//		}
	//		//再判断邻域索引子空间是否被当前层的支配，若是，加入聚类，若不是，不加入
	//		bool permit_overlap = true; //是否允许重叠聚类子空间
	//		std::vector<size_t> temp_layer;
	//		for (size_t i = 0; i < neigh_inx.size(); ++i) {
	//			bool add_in = false;
	//			if (permit_overlap|| flag_idx[neigh_inx[i]] == 0) {
	//				add_in = true;
	//			}
	//			else {
	//				add_in = false;
	//			}
	//			if (add_in) {
	//				size_t cur_rank = getMO_HLC().getSubspaceInfo(neigh_inx[i]).m_best_rank;
	//				std::vector<size_t> neigh_cur_compara;
	//				auto nei_inx = getMO_HLC().getSubspaceInfo(neigh_inx[i]).m_sub_neighbors;
	//				std::list<int> neighs;
	//				if (strict_nei) {
	//					//满足一定重合率的子空间邻域
	//					for (auto k : nei_inx) {
	//						auto center_box = getMO_HLC().subspaceTree().getBox(neigh_inx[i]);
	//						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//						size_t dim = INT16_MAX;
	//						for (size_t p = 0; p < center_box.size(); ++p) {
	//							auto bound1 = center_box[p];
	//							auto bound2 = compara_box[p];
	//							//判断该维是否重叠
	//							if (bound1.first < bound2.second && bound1.second > bound2.first) {
	//								dim = p;
	//								break;
	//							}
	//						}
	//						if (dim < INT16_MAX) {
	//							neighs.push_back(k);
	//						}
	//					}
	//				}
	//				else {
	//					for (auto k :nei_inx) {
	//						neighs.push_back(k);
	//					}
	//				}
	//				for (size_t j = 0; j < cur_cluster.size(); ++j) {
	//					if (std::find(neighs.begin(), neighs.end(), cur_cluster[j]) != neighs.end()) {
	//						neigh_cur_compara.push_back(cur_cluster[j]);
	//					}
	//				}
	//				for (size_t j = 0; j < neigh_cur_compara.size(); ++j) {//邻域子空间rank值比较
	//					size_t compara_rank = getMO_HLC().getSubspaceInfo(neigh_cur_compara[j]).m_best_rank;
	//					if (cur_rank >= compara_rank) {
	//						temp_layer.push_back(neigh_inx[i]);
	//						flag_idx[neigh_inx[i]] = 1;
	//						break;
	//					}
	//				}
	//			}
	//		}
	//		if (!temp_layer.empty()) {
	//			final_clusters[cluster_inx].emplace_back(temp_layer);
	//		}
	//		//最后更新当前类的最外层的rank
	//		int cur_layer_rank = -1*INT16_MAX;
	//		if (temp_layer.empty()) {
	//			cluster_cur_ranks[cluster_inx] = INT16_MAX;
	//		}
	//		else {
	//			for (size_t i = 0; i < temp_layer.size(); ++i) {
	//				int temp_rank = getMO_HLC().getSubspaceInfo(temp_layer[i]).m_best_rank;
	//				if (temp_rank > cur_layer_rank) {
	//					cur_layer_rank = temp_rank;
	//				}
	//			}
	//			cluster_cur_ranks[cluster_inx] = cur_layer_rank;
	//		}
	//	}
	//	//查看每个类的第一层是否被其他类包含，若是，则合并
	//	size_t num_eval = m_evaluations;
	//	m_his_effect_eval.push_back(num_eval);
	//	std::vector<bool> effective_flag(final_clusters.size(),true);
	//	for (size_t i = 0; i < final_clusters.size(); ++i) {
	//		auto front_layer = final_clusters[i][0][0];
	//		for (size_t j = 0; j < final_clusters.size(); ++j) {
	//			if (j != i) {
	//				bool flag=false;
	//				for (size_t k = 0; k < final_clusters[j].size(); ++k) {
	//					if (std::find(final_clusters[j][k].begin(), final_clusters[j][k].end(), front_layer) != final_clusters[j][k].end()) {
	//						effective_flag[i] = false;
	//						flag = true;
	//						break;
	//					}
	//				}
	//				if (flag) {
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	std::vector<std::vector<std::vector<size_t>>> final_effect_clusters;
	//	for (size_t i = 0; i < final_clusters.size(); ++i) {
	//		if (effective_flag[i]) {
	//			final_effect_clusters.emplace_back(final_clusters[i]);
	//		}
	//	}
	//	for (size_t i = 0; i < final_effect_clusters.size(); ++i) {
	//		getMO_HLC().getClusters().push_back(final_effect_clusters[i]);
	//	}
	//	//更新子空间所属类
	//	for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
	//		getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
	//	}
	//	for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
	//		for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
	//			for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
	//				getMO_HLC().getSubspaceInfo(getMO_HLC().getClusters()[i][j][k]).idx_cluster.push_back(i);
	//			}
	//		}
	//	}
	//}

	void SPMOEA1_0_1::clusterSubspace() {
		//每次聚类只加入一层，允许类间有重叠
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		std::vector<size_t> flag_idx(space_num, 0);
		std::vector<std::vector<std::vector<size_t>>> final_clusters;
		//先找出排序值最靠前的子空间
		while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
			//先得到现有子空间最好的排序值的索引
			int rank = INT16_MAX;//存放一个类最靠前的rank值
			std::vector<size_t> best_idx;//当前聚类具有最好rank值的子空间索引
			for (size_t i = 0; i < space_num; ++i) {
				if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank < rank) {
					rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
				}
			}
			for (size_t i = 0; i < space_num; ++i) {
				if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank == rank)
					best_idx.push_back(i);
			}
			//对这些最靠前的子空间进行聚类
			std::vector<std::vector<size_t>> cluster_idx;
			cluster_idx = clusterFrontSpace(best_idx);

			std::vector<std::vector<std::vector<size_t>>> cur_clusters;
			for (size_t i = 0; i < cluster_idx.size(); ++i) {
				std::vector<std::vector<size_t>> temp_clusters;
				temp_clusters.emplace_back(cluster_idx[i]);
				cur_clusters.emplace_back(temp_clusters);
			}

			for (size_t i = 0; i < cluster_idx.size(); ++i) {
				for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
					flag_idx[cluster_idx[i][j]] = 1;
				}
			}

			bool strict_neighbors = false;//是否重叠邻域
			bool strict_cluster = false;//对于每个子空间，只有其未聚类的邻域排序值全部小于该子空间，才加入到聚类
			bool permit_overlap = true;//是否允许聚类的子空间重叠
			auto clustered = cluster_idx;
			size_t add_count = cluster_idx.size();
			std::vector<bool> cluster_flag(cur_clusters.size(), true);
			while (add_count > 0) {
				add_count = 0;
				for (size_t i = 0; i < cur_clusters.size(); ++i) {
					std::vector<size_t> temp1;//每一个类的最外层
					if (cluster_flag[i]) {
						for (size_t j = 0; j < cur_clusters[i].back().size(); ++j) {
							size_t temp_rank = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_best_rank;
							auto temp_neighbor = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_sub_neighbors;
							std::vector<size_t> neighbors;
							for (auto& ss : temp_neighbor) {
								neighbors.push_back(ss);
							}
							if (strict_neighbors) {
								//满足一定重合率的子空间邻域
								std::vector<size_t> overlap_neighbor;
								for (auto& k : temp_neighbor) {
									auto center_box = getMO_HLC().subspaceTree().getBox(cur_clusters[i].back()[j]);
									auto compara_box = getMO_HLC().subspaceTree().getBox(k);
									size_t dim = INT16_MAX;
									for (size_t p = 0; p < center_box.size(); ++p) {
										auto bound1 = center_box[p];
										auto bound2 = compara_box[p];
										//判断该维是否重叠
										if (bound1.first < bound2.second && bound1.second > bound2.first) {
											dim = p;
											break;
										}
									}
									if (dim < INT16_MAX) {
										overlap_neighbor.push_back(dim);
									}
								}
								neighbors = overlap_neighbor;
							}
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
							for (size_t p = 0; p < sub.size(); ++p) {
								if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank >= temp_rank) {
									if (std::find(clustered[i].begin(), clustered[i].end(), sub[p]) == clustered[i].end()) {
										clustered[i].push_back(sub[p]);
										temp1.push_back(sub[p]);
										flag_idx[sub[p]] = 1;
										add_count++;
									}
								}
							}
						}
					}
					if (temp1.empty()) {
						cluster_flag[i] = false;
					}
					else {
						cur_clusters[i].emplace_back(temp1);
					}
				}
			}
			for (size_t i = 0; i < cur_clusters.size(); ++i) {
				final_clusters.emplace_back(cur_clusters[i]);
			}

			//while (!current.empty()) {
			//	std::vector<std::vector<size_t>> temp_current;//所有类的当前比较的子空间
			//	for (size_t i = 0; i < current.size(); ++i) {
			//		std::vector<size_t> temp1;//每一个类的第一层邻域
			//		if (current[i].size() > 1 || current[i][0] != INT16_MAX) {
			//			for (size_t j = 0; j < current[i].size(); ++j) {
			//				size_t temp_rank = getMO_HLC().getSubspaceInfo(current[i][j]).m_best_rank;
			//				auto temp_neighbor = getMO_HLC().getSubspaceInfo(current[i][j]).m_sub_neighbors;
			//				std::vector<size_t> neighbors;
			//				for (auto &ss : temp_neighbor) {
			//					neighbors.push_back(ss);
			//				}
			//				if (strict_neighbors) {
			//					//满足一定重合率的子空间邻域
			//					std::vector<size_t> overlap_neighbor;
			//					for (auto& k : temp_neighbor) {
			//						auto center_box = getMO_HLC().subspaceTree().getBox(current[i][j]);
			//						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
			//						size_t dim = INT16_MAX;
			//						for (size_t p = 0; p < center_box.size(); ++p) {
			//							auto bound1 = center_box[p];
			//							auto bound2 = compara_box[p];
			//							//判断该维是否重叠
			//							if (bound1.first < bound2.second && bound1.second > bound2.first) {
			//								dim = p;
			//								break;
			//							}
			//						}
			//						if (dim < INT16_MAX) {
			//							overlap_neighbor.push_back(dim);
			//						}
			//					}
			//					neighbors=overlap_neighbor;
			//				}
			//				std::vector<size_t> sub;
			//				for (auto& k : neighbors) {
			//					if (permit_overlap) {//是否允许聚类的子空间重叠
			//						sub.push_back(k);
			//					}
			//					else {
			//						if (flag_idx[k] == 0) {
			//							sub.push_back(k);
			//						}
			//					}
			//				}
			//				if (strict_cluster) {
			//					std::vector<Real> temp_space;
			//					bool add = true;
			//					for (size_t p = 0; p < sub.size(); ++p) {
			//						if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
			//							if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank > temp_rank) {
			//								temp_space.push_back(sub[p]);
			//							}
			//							else {
			//								add = false;
			//								break;
			//							}
			//						}
			//					}
			//					if (add) {
			//						for (auto& k : temp_space) {
			//							cluster_idx[i].push_back(k);
			//							temp1.push_back(k);
			//							flag_idx[k] = 1;
			//						}
			//					}
			//				}
			//				else {
			//					for (size_t p = 0; p < sub.size(); ++p) {
			//						if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank >= temp_rank) {
			//							if (std::find(cluster_idx[i].begin(), cluster_idx[i].end(), sub[p]) == cluster_idx[i].end()) {
			//								cluster_idx[i].push_back(sub[p]);
			//								temp1.push_back(sub[p]);
			//								flag_idx[sub[p]] = 1;
			//							}
			//						}
			//					}
			//				}
			//			}
			//		}
			//		if (temp1.empty()) {
			//			temp1.push_back(INT16_MAX);
			//		}
			//		cur_clusters[i].emplace_back(temp1);
			//		temp_current.push_back(temp1);
			//	}
			//	bool over = false;
			//	for (size_t i = 0; i < temp_current.size(); ++i) {
			//		if (temp_current[i][0] != INT16_MAX) {
			//			break;
			//		}
			//		else if (i == temp_current.size() - 1) {
			//			over = true;
			//		}
			//	}
			//	if (over) {
			//		current.clear();
			//		for (size_t i = 0; i < cur_clusters.size(); ++i) {
			//			if (cur_clusters[i].back()[0]==INT16_MAX) {
			//				cur_clusters[i].pop_back();
			//			}
			//			final_clusters.emplace_back(cur_clusters[i]);
			//		}
			//	}
			//	else {
			//		current = temp_current;
			//		temp_current.clear();
			//	}
			//}
		}
		for (size_t i = 0; i < final_clusters.size(); ++i) {
			getMO_HLC().getClusters().emplace_back(final_clusters[i]);
		}
		//更新子空间所属类
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
		}
		for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					getMO_HLC().getSubspaceInfo(getMO_HLC().getClusters()[i][j][k]).idx_cluster.push_back(i);
				}
			}
		}
	}

	//进行最优子空间的聚类：简单邻域聚类，还是考虑前沿解的分布？
	std::vector<std::vector<size_t>> SPMOEA1_0_1::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_0_1::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA1_0_1::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_0_1::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateVarSpaceInfo(temp_pop, pro, rnd);
	}

	void SPMOEA1_0_1::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_0_1::findSplitDim(int inx, Problem* pro) {
		auto& his_sols = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> front_ind_ratio;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			Real middle_v = (space_bound[j].first + space_bound[j].second) / 2;
			size_t upper_count = 0, lower_count = 0;
			for (size_t k = 0; k < his_sols.size(); ++k) {
				if (his_sols[k]->variable()[j] >= middle_v) {
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
			front_ind_ratio.push_back(min_v / max_v);
		}
		int dim = std::distance(front_ind_ratio.begin(), std::min_element(front_ind_ratio.begin(), front_ind_ratio.end()));

		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		return dim;
	}

	void SPMOEA1_0_1::NDSort(Population<Solution<>>& pop) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i) {
			objs.emplace_back(&pop[i].objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i) {
			pop[i].setFitness(rank[i]);
		}
	}

	void SPMOEA1_0_1::recordMetrics(Problem* pro, Algorithm* alg) {
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

	void SPMOEA1_0_1::initiObjSpace(Problem* pro) {
		//for (size_t i = 0; i < m_pop->size(); ++i) {
		//	m_pop->getOffspring()[i] = m_pop->at(i);
		//	m_pop->getParent()[i] = m_pop->at(i);
		//	m_pop->getOffspring()[i + m_pop->size()] = m_pop->at(i);
		//}
		///****************************/
		///*      计算目标值范围      */
		///****************************/
		////updateObjRange(pro);
		//size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		//m_pop_range.resize(obj_num);
		//for (int i = 0; i < obj_num; ++i) {
		//	m_pop_range[i].first = 1.0e14;
		//	for (int j = 0; j < m_pop->size(); ++j) {
		//		if (m_pop->at(j).objective()[i] < m_pop_range[i].first) {
		//			m_pop_range[i].first = m_pop->at(j).objective()[i];
		//		}
		//	}
		//}
		//for (int i = 0; i < obj_num; ++i) {
		//	m_pop_range[i].second = -1 * 1.0e14;
		//	for (int j = 0; j < m_pop->size(); ++j) {
		//		if (m_pop_range[i].second < m_pop->at(j).objective()[i]) {
		//			m_pop_range[i].second = m_pop->at(j).objective()[i];
		//		}
		//	}
		//}
		///****************************/
		///*      更新子目标最优      */
		///****************************/
		//for (size_t i = 0; i < obj_num; ++i) {
		//	m_subobj_opt_sol.emplace_back(m_pop->at(i));
		//}
		//updateSubObjOpt(m_pop->getParent());
		///****************************/
		///*    目标空间划分初始化    */
		///****************************/
		//for (size_t i = 0; i < m_num_region_obj; ++i) {
		//	m_obj_region_info.emplace_back(new ObjRegionInfo);
		//}
		//updateObjRange(m_pop->getOffspring(), pro);
		//m_front_pop_range = m_pop_range;
		//std::vector<std::pair<Real, Real>> m_obj_boundary;
		//if (m_normalize) {
		//	for (size_t i = 0; i < obj_num; ++i) {
		//		m_obj_boundary.push_back(std::make_pair<>(0., 1.));
		//	}
		//}
		//else {
		//	m_obj_boundary = getFrontPopRange();
		//}
		//m_obj_tree->setInitBox(m_obj_boundary);
		//m_obj_tree->inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		//m_obj_tree->buildIndex();

		//updateObjSpaceInfo(m_pop->getParent(), pro, false);
		//std::vector<size_t> space_idx;
		//for (size_t i = 0; i < numSubspace(); ++i) {
		//	auto& spaceinfo = getObjRegionInfo(i);
		//	if (spaceinfo.m_obj_rank == 0) {
		//		space_idx.push_back(i);
		//	}
		//}
		//m_front_space_inx.push_back(space_idx);
	}
}
