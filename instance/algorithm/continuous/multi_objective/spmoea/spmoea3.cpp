#include "spmoea3.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA3::initialize_() {
		SPMOEA::initialize_();
		/*auto& v = *m_param;
		size_t num_var = CAST_CONOP(m_problem.get())->numberVariables();
		size_t num_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		m_pop_size = v.get<int>("population size");
		if (m_pop_size % 2)
			throw MyExcept("Population size of NSGAII should be even.");
		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / num_var;
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
		m_num_region_var = v.get<int>("number of subspaces");
		m_num_region_obj = v.get<int>("number of obj subspaces");

		m_mo_hlc.emplace_back(new MO_HLC(num_var, num_obj));
		m_multi_pops.reset(new MultiPopulation<SPMOEA_pop>());

		std::vector<std::pair<Real, Real>> bound;
		for (size_t i = 0; i < num_var; ++i) {
			bound.emplace_back(CAST_CONOP(m_problem.get())->domain().range(i).limit);
		}

		initiVarBox(bound,m_problem.get());*/

		/*auto pro_name = CAST_CONOP(m_problem.get())->name();
		if (pro_name == "MOP_TEST3") {
			auto& v = GET_PARAM(m_problem->idParam());
			size_t peak_num = v.get<int>("numPeak");
			m_visit_count = std::vector<size_t>(peak_num);
		}*/
	}

	void SPMOEA3::initiVarBox(std::vector<std::pair<Real, Real>>& box, Problem* pro) {
		/*************************************/
		/*      MO_HLC搜索空间划分初始化     */
		/*************************************/
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		//getMO_HLC(0).initialVarSpace(box, m_num_region_var);//2次幂个子空间
		m_pop_var_range = box;
	}

	void SPMOEA3::run_() {
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

	void SPMOEA3::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		entry.push_back(m_R1);
		entry.push_back(m_R2);
		entry.push_back(m_R3);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(m_IGD);
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA3::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < getPop().size(); ++i)
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[0].push_back(&getPop()[i][j].phenotype());
				}
			/*for (size_t i = 0; i < m_pop->getOffspring().size(); ++i)
				m_solution[1].push_back(&m_pop->getOffspring()[i].phenotype());*/
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA3::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//每个子空间初始化一个子种群
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
		size_t space_num = getMO_HLC().numSubspace();
		for (size_t i = 0; i < space_num; ++i) {
			SPMOEA_pop temp_pop(pop_num, pro);
			std::vector<std::vector<Real>> sols;
			auto bound = getMO_HLC().subspaceTree().getBox(i);
			for (size_t j = 0; j < pop_num; ++j) {
				for (size_t k = 0; k < bound.size(); ++k) {
					temp_pop[j].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
				}
				temp_pop[j].evaluate(pro, alg);
				sols.emplace_back(temp_pop[j].variable().vect());
			}
			NDSort(temp_pop);
			//初始化子代
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring()[j] = temp_pop[j];
				temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			}
			temp_pop.setRate(getCr(), getMr());
			temp_pop.setEta(getCeta(), getMeta());
			getPop().append(temp_pop);
			//更新种群信息
			getPop().back().setSearchBox(bound);
			//更新子空间信息
			getMO_HLC().getSubspaceInfo(i).m_feature_state = 1;
			updateVarSpaceInfo(temp_pop, pro, rnd);
			updateHistoryFrontSol(temp_pop, pro);
			SPMOEA::updateArchive(archiveNum(), pro);
			SPMOEA::updateNewPop(temp_pop);
			SPMOEA::updateHistorySols(temp_pop);
			SPMOEA::updateObjRange(temp_pop, pro);
			updateObjSpace();
			////根据子空间排序值聚类
			//clusterSubspace();

			SPMOEA::recordMetrics(pro, alg);
			//SPMOEA::record();
			//SPMOEA::recordMetrics(pro, alg);
		}
		updateEE(space_num * pop_num, 0);


		////m_pop.reset(new SPMOEA1_pop(m_pop_size, m_problem.get()));
		////initiVarSpace(m_problem.get());
		////在子空间生成子种群
		//for (size_t i = 0; i < getMO_HLC(0).numSubspace(); ++i) {
		//	SPMOEA_pop temp_pop(m_sub_pop_size, pro);
		//	//设置初始化区域
		//	auto sub_bound = getMO_HLC(0).subspaceTree().getBox(i);
		//	CAST_CONOP(pro)->setInitialDomain(sub_bound);//不能修改问题初始边界
		//	temp_pop.initialize(pro,rnd);
		//	temp_pop.evaluate(pro, alg);
		//	NDSort(temp_pop);
		//	temp_pop.updatePopDistribute(pro);
		//	temp_pop.addPopdist();
		//	getPop().append(temp_pop);
		//	
		//	
		//	getPop()[i].setRate(getCr(), getMr());
		//	getPop()[i].setEta(getCeta(), getMeta());
		//	getPop()[i].setSearchRange(sub_bound);
		//	////输出位置
		//	//std::cout << "第" << i << "个种群位置" << std::endl;
		//	//for (size_t j = 0; j < temp_pop.size(); ++j) {
		//	//	std::cout << temp_pop[j].variable().vect()[0] << "   " << temp_pop[j].variable().vect()[1] << std::endl;
		//	//}
		//}
		////int tag = EvaluationTag::kNormalEval;
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	//初始化子种群子代
		//	for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//		getPop()[i].getOffspring()[j] = getPop()[i][j];
		//		getPop()[i].getOffspring()[j + getPop()[i].size()] = getPop()[i][j];
		//	}
		//}
		//SPMOEA_pop temp_pop(m_sub_pop_size* getPop().size(), pro);
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	for (size_t j= 0; j < getPop()[i].size(); ++j) {
		//		temp_pop[i*m_sub_pop_size+j] = getPop()[i][j];
		//	}
		//}
		//NDSort(temp_pop);
		////更新子空间信息
		//SPMOEA::updateVarSpaceInfo(temp_pop, pro, rnd);
		////更新历史前沿
		//SPMOEA::updateHistoryFrontSol(temp_pop, pro);
		//SPMOEA::updateArchive(archiveNum(), pro);
		//SPMOEA::updateObjRange(temp_pop,pro);
		//updateObjSpace();
		//
		//initiObjSpace(m_problem.get());
	}

	void SPMOEA3::initiObjSpace(Problem* pro) {
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

	int SPMOEA3::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		//根据历史解分布判断是否存在局部结构
		for (size_t i = 0; i < getPop().size(); ++i) {
			auto search_box = getPop()[i].getPopSearchBox();

		}


		//若存在，划分出局部结构，重新布置种群搜索

		//若不存在，根据历史解的分布制定采样策略

		//更新搜索空间历史解、存档，记录指标

		//更新子种群



		//检查各子种群收敛位置，如果判定为收敛到边缘，将前沿所在区域划分出来，标记“exploit”子空间
		std::vector<bool> pop_converge_flag;
		for (size_t i = 0; i < getPop().size(); ++i) {
			pop_converge_flag.push_back(popConverged(i));
		}
		for (size_t i = 0; i < pop_converge_flag.size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				if (pop_converge_flag[i]) {
					auto indicator = subspaceSeperable(i);
					if (!indicator) {
						//或者每次只对空间最大的子空间细分
						//将该子空间种群收敛范围区分出来，重新划分
						//先保存子空间历史信息
						if (m_cur_split_times < m_max_split_times) {//还是使用容积限制条件禁止无限分解？
							auto hist_info = getMO_HLC().getSubspaceInfo(i).m_history_inds;
							//再细分子空间
							size_t pre_subspace_num = getMO_HLC().numSubspace();
							splitSubspace(i, 1);
							m_cur_split_times++;
							//重新分配第i个子种群，
							auto bound = getMO_HLC().subspaceTree().getBox(i);
							updatePop(i, bound, pro, alg, rnd);
							//在新的子空间生成子种群
							std::vector<size_t> add_space_inx;
							for (size_t j = pre_subspace_num; j < getMO_HLC().numSubspace(); ++j) {
								add_space_inx.push_back(j);
							}
							size_t pop_num = getPop().size();
							for (size_t j = 0; j < add_space_inx.size(); ++j) {
								auto bound = getMO_HLC().subspaceTree().getBox(add_space_inx[j]);
								generatePop(bound, pro, alg, rnd);
							}
							//使用新产生的解更新子空间信息
							SPMOEA_pop temp_pop(m_sub_pop_size * (add_space_inx.size() + 1), pro);
							for (size_t j = 0; j < m_sub_pop_size; ++j) {
								temp_pop[j] = getPop()[i][j];
							}
							for (size_t j = pop_num; j < getPop().size(); ++j) {
								for (size_t k = 0; k < getPop()[j].size(); ++k) {
									temp_pop[(j - pop_num + 1) * m_sub_pop_size + k] = getPop()[j][k];
								}
							}
							NDSort(temp_pop);
							SPMOEA::updateVarSpaceInfo(temp_pop, pro, rnd);
							SPMOEA::updateHistoryFrontSol(temp_pop, pro);
							SPMOEA::updateArchive(archiveNum(), pro);
							SPMOEA::updateObjRange(temp_pop, pro);
							//更新第i个子空间历史信息
							SPMOEA::updateHisVarSpaceInfo(hist_info, pro, rnd);
						}
					}
					else {
						//对该子空间对半分，
						getPop()[i].setPopState("converged");
					}
				}
			}
		}

		//子种群演化
		bool m_evolve_by_prediction = false;
		if (!m_evolve_by_prediction) {
			//generate solutions by operators, such as GA, DE
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (getPop()[i].getPopState() == "active") {
					generateOffspring(pro, rnd);
				}
			}
		}
		else {
			//generating solutions by learning methods
			bool m_prediction = false;
			if (!m_prediction) {
				//multi-population method and operators
				//子空间预采样，更新子空间前沿解，子空间内代表个体排序聚类，多种群搜索

			}
			else {
				//multi-popution and prediction
				//代内排序预测，代际预测更新
				//子种群演化历史轨迹
			}
		}

		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				SPMOEA_pop offspring_pop(getPop()[i].size(), pro);
				for (size_t j = 0; j < getPop()[i].size(); j++) {
					//评价子代
					tag = getPop()[i].getOffspring()[j].evaluate(pro, alg);
					if (tag != EvaluationTag::kNormalEval)
						break;
					offspring_pop[j] = getPop()[i].getOffspring()[j];
				}
				//更新子空间信息
				NDSort(offspring_pop);
				updateVarSpaceInfo(offspring_pop, pro, rnd);
				//更新历史信息
				updateHistoryFrontSol(offspring_pop, pro);
			}
		}

		//子种群内部排序、淘汰选择
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
				}
				//m_multi_pop->at(i).envirSelection(pro, rnd, 0);
				getPop()[i].envirSelection(getPop()[i], getPop()[i].getOffspring());

				////输出位置
				//std::cout << "第" << i << "个种群位置" << std::endl;
				//for (size_t j = 0; j < getPop()[i].size(); ++j) {
				//	std::cout << getPop()[i][j].variable().vect()[0] << "   " << getPop()[i][j].variable().vect()[1] << std::endl;
				//}

				updatePopDistribute(i, pro);
				updatePopdist(i);
			}
		}

		//根据各个子种群或子空间的搜索结果更新外部存档,
		//使用历史所有非支配解更新archive
		updateArchive(archiveNum(), pro);

		recordMetrics(pro, alg);

		return tag;
	}

	void SPMOEA3::generatePop(std::vector<std::pair<Real, Real>>& bound, Problem* pro, Algorithm* alg, Random* rnd) {
		SPMOEA_pop temp_pop(m_sub_pop_size, pro);
		CAST_CONOP(pro)->setInitialDomain(bound);
		temp_pop.initialize(pro, rnd);
		temp_pop.evaluate(pro, alg);
		//NDSort(temp_pop);
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop.getOffspring()[i] = temp_pop[i];
			temp_pop.getOffspring()[i + temp_pop.size()] = temp_pop[i];
		}
		//更新种群状态
		temp_pop.updatePopDistribute(pro);
		temp_pop.addPopdist();

		temp_pop.setRate(getCr(), getMr());
		temp_pop.setEta(getCeta(), getMeta());
		temp_pop.setSearchRange(bound);

		getPop().append(temp_pop);

	}

	void SPMOEA3::updatePop(size_t inx, std::vector<std::pair<Real, Real>>& bound, Problem* pro, Algorithm* alg, Random* rnd) {
		SPMOEA_pop temp_pop(m_sub_pop_size, pro);
		CAST_CONOP(pro)->setInitialDomain(bound);
		temp_pop.initialize(pro, rnd);
		temp_pop.evaluate(pro, alg);
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			getPop()[inx][i] = temp_pop[i];
		}
		//更新种群状态
		getPop()[inx].updatePopDistribute(pro);
		getPop()[inx].addPopdist();

		getPop()[inx].setRate(getCr(), getMr());
		getPop()[inx].setEta(getCeta(), getMeta());
		getPop()[inx].setSearchRange(bound);

	}

	void SPMOEA3::generateOffspring(Problem* pro, Random* rnd) {
		auto& m_pop = getPop();
		for (size_t i = 0; i < m_pop.size(); ++i) {
			auto sub_bound = getMO_HLC().subspaceTree().getBox(i);
			for (size_t j = 0; j < m_pop[i].size(); j += 2) {
				std::vector<size_t> p(2);
				do {//不重复采样
					do {
						p[0] = m_pop[i].tournamentSelection(pro, rnd);
						p[1] = m_pop[i].tournamentSelection(pro, rnd);
					} while (p[1] == p[0]);
					m_pop[i].crossover(p[0], p[1], m_pop[i].getOffspring()[j], m_pop[i].getOffspring()[j + 1], pro, rnd);
					m_pop[i].mutate(m_pop[i].getOffspring()[j], pro, rnd);
					m_pop[i].mutate(m_pop[i].getOffspring()[j + 1], pro, rnd);
					//子代越界处理
					repairSol(m_pop[i].getOffspring()[j].variable().vect(), sub_bound, rnd);
					repairSol(m_pop[i].getOffspring()[j + 1].variable().vect(), sub_bound, rnd);
				} while (m_pop[i].getOffspring()[j].variable().vect() == m_pop[i][p[0]].variable().vect() || \
					m_pop[i].getOffspring()[j].variable().vect() == m_pop[i][p[1]].variable().vect() || \
					m_pop[i].getOffspring()[j + 1].variable().vect() == m_pop[i][p[0]].variable().vect() || \
					m_pop[i].getOffspring()[j + 1].variable().vect() == m_pop[i][p[1]].variable().vect());

				m_pop[i].getOffspring()[j].setCounter(0);
				m_pop[i].getOffspring()[j + 1].setCounter(0);
			}
			////输出种群
			//for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
			//	std::cout << m_multi_pop->at(i).getOffspring()[j].variable()[0] << "      " << m_multi_pop->at(i).getOffspring()[j].variable()[1] << std::endl;
			//}
			//std::cout << "*******************" << std::endl;
		}
	}

	bool SPMOEA3::updateCluster() {
		return false;
	}

	void SPMOEA3::clusterSubspace() {
		getMO_HLC().rankClustering();
	}

	void SPMOEA3::updateHistoryFrontSol(Population<Solution<>>& pop, Problem* pro) {
		//先加入每代的前沿个体
		//pop的排序值不对，找出pop中的非支配解集
		//NDSort(pop);
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0) {
				temp_pop.append(pop[i]);
			}
		}
		m_gen_front_sols.emplace_back(std::make_shared<Population<Solution<>>>(temp_pop));
		//再更新历史非支配解
		if (m_history_front_sols.empty()) {
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				m_history_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
			}
		}
		else {
			std::vector<size_t> temp_his_sols(m_history_front_sols.size(), 1);
			std::vector<size_t> temp_pop_sols(temp_pop.size(), 0);
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
					if (temp_his_sols[i] == 1) {
						Dominance dominanceship = temp_pop[j].compare(*m_history_front_sols[i], CAST_CONOP(pro)->optimizeMode());
						if (dominanceship == Dominance::kNonDominated) {
							if (i == m_history_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else if (dominanceship == Dominance::kDominant) {
							temp_his_sols[i] = 0;
							if (i == m_history_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else {
							break;
						}
					}
				}
			}
			auto temp = m_history_front_sols;
			m_history_front_sols.clear();
			for (size_t i = 0; i < temp.size(); ++i) {
				if (temp_his_sols[i] == 1) {
					m_history_front_sols.emplace_back(temp[i]);
				}
			}
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				if (temp_pop_sols[i] == 1) {
					m_history_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
				}
			}
		}
	}

	void SPMOEA3::updateArchive(size_t num, Problem* pro) {
		//从历史前沿中选择一定数量的个体作为archive
		//三种选择方式：1、拥挤选择；2、参考向量选择；3、子空间跨度选择
		if (m_history_front_sols.size() <= 2 * num) {
			m_archive.clear();
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				m_archive.emplace_back(m_history_front_sols[i]);
			}
		}
		else {
			//目标子空间选择
			std::vector<size_t> select_index = selectIndiFromFront(m_history_front_sols, num, pro);
			m_archive.clear();
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				if (std::find(select_index.begin(), select_index.end(), i) != select_index.end()) {
					m_archive.emplace_back(m_history_front_sols[i]);
				}
			}
			std::cout << "累积前沿解个数" << m_history_front_sols.size() << std::endl;
		}
	}

	int SPMOEA3::subspaceSeperable(size_t inx) {
		//计算到子空间边缘的最小距离平均距离
		auto bound = getPop()[inx].getSearchRange();
		std::vector<Real> dim_span;
		for (size_t i = 0; i < bound.size(); ++i) {
			dim_span.push_back(bound[i].second - bound[i].first);
		}
		auto min_span = *std::min_element(dim_span.begin(), dim_span.end());

		auto min_bound_dist = minDistFromBound(inx);
		Real mean_dist = 0.;
		for (size_t i = 0; i < min_bound_dist.size(); ++i) {
			mean_dist += min_bound_dist[i];
		}
		mean_dist /= min_bound_dist.size();
		if (mean_dist / min_span < 0.1) {
			return 1;
		}
		else {
			return 0;
		}
	}

	std::vector<std::pair<Real, Real>> SPMOEA3::distFromBound(size_t pop_inx) {
		//点到边界的最小距离和
		auto bound = getPop()[pop_inx].getSearchRange();
		std::vector<std::pair<Real, Real>> dist;
		for (size_t i = 0; i < bound.size(); ++i) {
			Real sum_d = 0;
			for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
				sum_d += getPop()[pop_inx][j].variable().vect()[i];
			}
			auto mean_d = sum_d / getPop()[pop_inx].size();
			dist.emplace_back(std::make_pair<>(std::fabs(mean_d - bound[i].first), std::fabs(mean_d - bound[i].second)));
		}
		return dist;
	}

	std::vector<Real> SPMOEA3::minDistFromBound(size_t pop_inx) {
		//点到边界的最小距离和
		auto bound = getPop()[pop_inx].getSearchRange();
		std::vector<Real> min_dist;
		for (size_t i = 0; i < getPop()[pop_inx].size(); ++i) {
			auto temp = (Real)INT16_MAX;
			auto var = getPop()[pop_inx][i].variable().vect();
			for (size_t j = 0; j < bound.size(); ++j) {
				if (temp > std::fabs(bound[j].first - var[j])) {
					temp = std::fabs(bound[j].first - var[j]);
				}
				if (temp > std::fabs(bound[j].second - var[j])) {
					temp = std::fabs(bound[j].second - var[j]);
				}
			}
			min_dist.push_back(temp);
		}
		return min_dist;
	}

	void SPMOEA3::splitSubspace(size_t inx, int flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}
}
