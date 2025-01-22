#include "SPMOEA1_4.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_4::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_4::run_() {
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

	void SPMOEA1_4::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_4::updateBuffer() {
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

	void SPMOEA1_4::initiVarSpace(Problem *pro) {
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

	void SPMOEA1_4::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
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
		temp_pop.initialize(pro,rnd);
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
			getMO_HLC().updateSubspaceInfo(temp_off_pop, pro, rnd);
		}
		updateVarSpace(pro, rnd);
		updateEE(pop_num, 0);
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
	}

	int SPMOEA1_4::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		/********************************************************************************
							   根据子种群前沿子空间的体积分配搜索资源
		********************************************************************************/
		//单纯使用子种群信息属于不完全信息比较，可否使用子空间全部的有效信息进行子空间的比较？
		//对所有子种群个体进行排序
		//根据前面的front_in_bound_ratio确定每个子种群的搜索资源和探索开发比
		std::vector<size_t> assign_pop_resource;
		size_t eval_num = 200;
		PopResourceAssign(assign_pop_resource, eval_num, pro);
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
			for (size_t i = 0; i < getPop().size(); ++i) {
				/*if (m_evaluations < eval_num) {
					interactive_type.push_back(1);
				}
				else {
					interactive_type.push_back(2);
				}*/
				if (m_divide_iteration % 2 == 0) {
					interactive_type.push_back(1);
				}
				else {
					interactive_type.push_back(2);
				}
			}
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
			Population<Solution<>> new_pop;
			for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); j++) {
				new_pop.append(getPop()[i].getOffspring()[j]);
			}
			getPop()[i].updatePopHisInfo(new_pop, pro);
			//getMO_HLC().updateSubspaceInfo(new_pop, pro, rnd);//子代更新子空间前沿解和代表解
		}

		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		if (m_divide_iteration % 3 == 0) {
			getMO_HLC().updateSubspaceInfo(offspring_pop, pro, rnd);//子代更新子空间前沿解和代表解
			updateVarSpace(pro, rnd);//使用子种群子代更新子空间信息
			spaceSubdivision(pro, rnd);
			updateFrontSpace();
			clusterSubspace();
		}

		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		SPMOEA::NDSort(offspring_pop);
		
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//使用历史所有非支配解更新archive
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space= getMO_HLC().numSubspace();
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				auto v = getMO_HLC().subspaceTree().getBoxVolume(i);
				front_ratio += (v / total_volume);
				front_count++;
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
			auto& bound = getPop()[i].getPopSearchBox();
			getPop()[i].envirSelection(pro, rnd, 0);
			getPop()[i].iteration()++;

			updatePopDistribute(i, pro);
			updatePopdist(i);
		}
		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA1_4::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
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

	bool SPMOEA1_4::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_4::updateSubPopSpace(size_t pop_inx, Problem *pro, Random *rnd) {
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

	void SPMOEA1_4::updateVarSpace(Problem *pro, Random *rnd) {
		//更新子空间的最好排序值，取子空间的代表个体还是前沿个体？
		std::vector<int> pre_rank;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			pre_rank.push_back(getMO_HLC().getSubspaceInfo(i).m_best_rank);
		}
		std::vector<std::shared_ptr<Solution<>>> temp_pop;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getSubspaceInfo(i).m_represent_sol.size(); ++j) {
				temp_pop.emplace_back(getMO_HLC().getSubspaceInfo(i).m_represent_sol[j]);
			}
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			objs.emplace_back(&temp_pop[i]->objective());
		}
		std::vector<int> rank;
		int num_layer = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		getMO_HLC().setLayerNum(num_layer);
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop[i]->setFitness(rank[i]);
		}
		//更新子空间最好排序值
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).m_best_rank = INT16_MAX;
		}
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			auto temp_sol = temp_pop[i]->variable().vect();
			auto inx = getMO_HLC().subspaceTree().getRegionIdx(temp_sol);
			if (temp_pop[i]->fitness() < getMO_HLC().getSubspaceInfo(inx).m_best_rank) {
				getMO_HLC().getSubspaceInfo(inx).m_best_rank = temp_pop[i]->fitness();
			}
		}
		//对于没有个体的子空间，找邻域最好的排序值
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == INT16_MAX) {
				auto neigh_space = getMO_HLC().getSubspaceInfo(i).m_sub_neighbors;
				std::vector<int> neigh_ranks;
				auto iter = neigh_space.begin();
				for (size_t j = 0; j < neigh_space.size(); ++j) {
					if (getMO_HLC().getSubspaceInfo(*iter).m_best_rank < INT16_MAX) {
						neigh_ranks.push_back(getMO_HLC().getSubspaceInfo(*iter).m_best_rank);
					}
					++iter;
				}
				if (neigh_ranks.empty()) {
					getMO_HLC().getSubspaceInfo(i).m_best_rank = num_layer;
				}
				else if (neigh_ranks.size() == 1) {
					getMO_HLC().getSubspaceInfo(i).m_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
				}
				else {
					int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
					int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
					getMO_HLC().getSubspaceInfo(i).m_best_rank = std::ceil((min_rank + max_rank) / 2) + 1;
				}
			}
		}

		std::vector<int> past_rank;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			past_rank.push_back(getMO_HLC().getSubspaceInfo(i).m_best_rank);
		}
	}

	void SPMOEA1_4::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t eval_num, Problem *pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		updateFrontRegionLinkSpace();
		if (m_divide_iteration % 2 == 0) {
			for (size_t i = 0; i < getPop().size(); ++i) {
				assign_pop_resource.push_back(max_pop_size);
			}
		}
		else {
			auto front_link_spaces = getFrontRegionLinkSpace();
			if (front_link_spaces.size() == 1) {
				assign_pop_resource.push_back(max_pop_size);
			}
			else {
				assign_pop_resource.push_back(7 * front_link_spaces.size());
			}
		}
		/*if (m_evaluations < eval_num) {
			for (size_t i = 0; i < getPop().size(); ++i) {
				assign_pop_resource.push_back(max_pop_size);
			}
		}
		else {
			updateFrontRegionLinkSpace();
			auto front_link_spaces = getFrontRegionLinkSpace();
			if (front_link_spaces.size() == 1) {
				assign_pop_resource.push_back(max_pop_size);
			}
			else {
				assign_pop_resource.push_back(3 * front_link_spaces.size());
			}
		}*/
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_4::generateOffspring(Problem *pro, Random *rnd, const std::vector<size_t>& pop_resource, std::vector<int> type) {
		size_t num_exploit = 0;
		auto search_bound = CAST_CONOP(pro)->boundary();
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
			//在子区域前排个体间开发,开发方式选择与邻域个体之间交互
			if (type[i] == 1) {
				//auto all_off = sampleByGA(getPop()[i], exploit_num, pro, rnd);
				size_t kk = 2;
				auto bound = CAST_CONOP(pro)->boundary();
				auto all_off = sampleByDE(getPop()[i],bound, exploit_num, kk,pro, rnd);
				//为子代赋值
				for (size_t k = 0; k < all_off.size(); ++k) {
					getPop()[i].getOffspring()[k].variable().vect() = all_off[k];
					getPop()[i].getOffspring()[k].setCounter(0);
				}
			}
			else if (type[i] == 2) {//前沿子空间连通集交互，连通集中采样
				auto front_link_spaces = getFrontRegionLinkSpace();
				if (front_link_spaces.size() == 1) {//子空间连通
					auto manifold_dist = calSpaceManifoldDist(front_link_spaces[0]);//优化计算
					//两端扩展,分别采样一个解
					auto off1 = sampleByLinkSpaces(front_link_spaces[0], manifold_dist,2, pro, rnd);
					//子代赋值
					getPop().back().getOffspring()[0].variable().vect() = off1[0];
					getPop().back().getOffspring()[1].variable().vect() = off1[1];
					getPop().back().getOffspring()[0].setCounter(0);
					getPop().back().getOffspring()[1].setCounter(0);

					/*中间采样*/
					//连通子空间内采样：1、随机选择一个子空间，与其邻域子空间交互；2、根据解在目标空间的分布交互
					//先获取连通子空间内的前沿解在目标空间的分布，然后根据目标空间前沿聚类之后的稀疏度采样
					auto all_off = sampleByManifoldSpaces(front_link_spaces[0], manifold_dist, pop_resource.back() - 2, pro, rnd);
					for (size_t j = 0; j < all_off.size(); ++j) {
						getPop().back().getOffspring()[2 + j].variable().vect() = all_off[j];
						getPop().back().getOffspring()[2 + j].setCounter(0);
					}
				}
				else {
					//使子连通集扩展
					std::vector<std::vector<Real>> all_off;
					size_t off_count = 0;
					for (size_t j = 0; j < front_link_spaces.size(); ++j) {
						auto manifold_dist = calSpaceManifoldDist(front_link_spaces[j]);
						auto off = sampleByLinkSpaces(front_link_spaces[j], manifold_dist, 2,pro, rnd);
						//为子代赋值
						for (size_t k = 0; k < off.size(); ++k) {
							getPop().back().getOffspring()[off_count].variable().vect() = off[k];
							getPop().back().getOffspring()[off_count].setCounter(0);
							off_count++;
						}
					}
					//子连通集之间的交互（填补空缺）
					auto select_pair_spaces = findCloseSpaces(front_link_spaces, rnd);
					//分别在每一个子连通集上进行最小距离子空间的交互
					for (size_t j = 0; j < select_pair_spaces.size(); ++j) {
						size_t inx1 = std::get<0>(select_pair_spaces[j]);
						size_t inx2 = std::get<1>(select_pair_spaces[j]);
						//在这两个子空间采用向量方式采样
						//子空间前沿个体交互
						auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
						auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
						size_t s1 = std::floor(front_sols1.size() * rnd->uniform.next());
						size_t s2 = std::floor(front_sols2.size() * rnd->uniform.next());
						auto p1 = front_sols1[s1]->variable().vect();
						auto p2 = front_sols2[s2]->variable().vect();
						////子空间随机选择交互
						//std::vector<Real> p1;
						//std::vector<Real> p2;
						//auto& box1 = getMO_HLC().subspaceTree().getBox(inx1);
						//auto& box2 = getMO_HLC().subspaceTree().getBox(inx2);
						//for (size_t k = 0; k < box1.size(); ++k) {
						//	p1.push_back(box1[k].first + (box1[k].second - box1[k].first) * rnd->uniform.next());
						//	p2.push_back(box2[k].first + (box2[k].second - box2[k].first) * rnd->uniform.next());
						//}
						
						for (size_t k = 0; k < 5; ++k) {
							auto off = vectorOperator(p1, p2, search_bound,0.1, rnd);
							getPop().back().getOffspring()[off_count].variable().vect() = off;
							getPop().back().getOffspring()[off_count].setCounter(0);
							off_count++;
						}
					}
				}
			}
		}
		updateEE(0, num_exploit);
	}

	void SPMOEA1_4::spaceSubdivision(Problem *pro, Random *rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		std::vector<size_t> front_spaces;
		for (size_t i = 0; i < num_spaces; ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_spaces.push_back(i);
			}
		}
		//细分子空间
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			size_t min_num = getSplitGranularity();
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

		//子空间体积占比
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_4::clusterSubspace() {
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

		///*********前沿子空间逐层聚**********/
		////每次聚类只加入一层，允许类间有重叠
		//getMO_HLC().getClusters().clear();
		//auto space_num = getMO_HLC().numSubspace();
		//std::vector<size_t> flag_idx(space_num, 0);
		//std::vector<std::vector<std::vector<size_t>>> final_clusters;
		////先找出排序值最靠前的子空间
		//while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
		//	//先得到现有子空间最好的排序值的索引
		//	int rank = INT16_MAX;//存放一个类最靠前的rank值
		//	std::vector<size_t> best_idx;//当前聚类具有最好rank值的子空间索引
		//	for (size_t i = 0; i < space_num; ++i) {
		//		if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank < rank) {
		//			rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
		//		}
		//	}
		//	for (size_t i = 0; i < space_num; ++i) {
		//		if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank == rank)
		//			best_idx.push_back(i);
		//	}
		//	//对这些最靠前的子空间进行聚类
		//	std::vector<std::vector<size_t>> cluster_idx;
		//	cluster_idx = clusterFrontSpace(best_idx);

		//	std::vector<std::vector<std::vector<size_t>>> cur_clusters;
		//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
		//		std::vector<std::vector<size_t>> temp_clusters;
		//		temp_clusters.emplace_back(cluster_idx[i]);
		//		cur_clusters.emplace_back(temp_clusters);
		//	}

		//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
		//		for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
		//			flag_idx[cluster_idx[i][j]] = 1;
		//		}
		//	}

		//	bool strict_neighbors = false;//是否重叠邻域
		//	bool strict_cluster = false;//对于每个子空间，只有其未聚类的邻域排序值全部小于该子空间，才加入到聚类
		//	bool permit_overlap = true;//是否允许聚类的子空间重叠
		//	auto clustered = cluster_idx;
		//	size_t add_count = cluster_idx.size();
		//	std::vector<bool> cluster_flag(cur_clusters.size(), true);
		//	while (add_count > 0) {
		//		add_count = 0;
		//		for (size_t i = 0; i < cur_clusters.size(); ++i) {
		//			std::vector<size_t> temp1;//每一个类的最外层
		//			if (cluster_flag[i]) {
		//				for (size_t j = 0; j < cur_clusters[i].back().size(); ++j) {
		//					size_t temp_rank = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_best_rank;
		//					auto temp_neighbor = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_sub_neighbors;
		//					std::vector<size_t> neighbors;
		//					for (auto& ss : temp_neighbor) {
		//						neighbors.push_back(ss);
		//					}
		//					if (strict_neighbors) {
		//						//满足一定重合率的子空间邻域
		//						std::vector<size_t> overlap_neighbor;
		//						for (auto& k : temp_neighbor) {
		//							auto center_box = getMO_HLC().subspaceTree().getBox(cur_clusters[i].back()[j]);
		//							auto compara_box = getMO_HLC().subspaceTree().getBox(k);
		//							size_t dim = INT16_MAX;
		//							for (size_t p = 0; p < center_box.size(); ++p) {
		//								auto bound1 = center_box[p];
		//								auto bound2 = compara_box[p];
		//								//判断该维是否重叠
		//								if (bound1.first < bound2.second && bound1.second > bound2.first) {
		//									dim = p;
		//									break;
		//								}
		//							}
		//							if (dim < INT16_MAX) {
		//								overlap_neighbor.push_back(dim);
		//							}
		//						}
		//						neighbors = overlap_neighbor;
		//					}
		//					std::vector<size_t> sub;
		//					for (auto& k : neighbors) {
		//						if (permit_overlap) {//是否允许聚类的子空间重叠
		//							sub.push_back(k);
		//						}
		//						else {
		//							if (flag_idx[k] == 0) {
		//								sub.push_back(k);
		//							}
		//						}
		//					}
		//					for (size_t p = 0; p < sub.size(); ++p) {
		//						if (getMO_HLC().getSubspaceInfo(sub[p]).m_best_rank >= temp_rank) {
		//							if (std::find(clustered[i].begin(), clustered[i].end(), sub[p]) == clustered[i].end()) {
		//								clustered[i].push_back(sub[p]);
		//								temp1.push_back(sub[p]);
		//								flag_idx[sub[p]] = 1;
		//								add_count++;
		//							}
		//						}
		//					}
		//				}
		//			}
		//			if (temp1.empty()) {
		//				cluster_flag[i] = false;
		//			}
		//			else {
		//				cur_clusters[i].emplace_back(temp1);
		//			}
		//		}
		//	}
		//	for (size_t i = 0; i < cur_clusters.size(); ++i) {
		//		final_clusters.emplace_back(cur_clusters[i]);
		//	}
		//}
		//for (size_t i = 0; i < final_clusters.size(); ++i) {
		//	getMO_HLC().getClusters().emplace_back(final_clusters[i]);
		//}
		////更新子空间所属类
		//for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//	getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
		//}
		//for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
		//	for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
		//		for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
		//			getMO_HLC().getSubspaceInfo(getMO_HLC().getClusters()[i][j][k]).idx_cluster.push_back(i);
		//		}
		//	}
		//}
	}

	std::vector<std::vector<size_t>> SPMOEA1_4::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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

	bool SPMOEA1_4::subspaceLink(size_t inx1, size_t inx2) {
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

	void SPMOEA1_4::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_4::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem *pro, Random *rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateVarSpaceInfo(temp_pop, pro, rnd);
	}

	void SPMOEA1_4::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_4::findSplitDim(int inx, Problem *pro) {
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

	void SPMOEA1_4::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_4::recordMetrics(Problem *pro, Algorithm *alg) {
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

	void SPMOEA1_4::initiObjSpace(Problem *pro) {

	}
}