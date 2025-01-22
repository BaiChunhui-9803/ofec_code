#include "spmoea9.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/functional.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include<fstream>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA9::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_pop_size = v.get<int>("population size");
		if (m_pop_size % 2)
			throw MyExcept("Population size of NSGAII should be even.");
		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
		m_num_region_var = v.get<int>("number of subspaces");
		m_num_region_obj = v.get<int>("number of obj subspaces");

		m_add_neighbor = v.get<bool>("cluster add neighbors");
		m_sample_in_basin = v.get<bool>("sample in basin");
		m_subspace_select = v.get<bool>("select by subspace");
		m_evolve_by_potential = v.get<bool>("evolve by potential");
		m_mean_potential = v.get<bool>("mean basin potential");

		m_mo_hlc.reset(new MO_HLC(CAST_CONOP(m_problem.get())->numberVariables(), CAST_CONOP(m_problem.get())->numberObjectives()));
		//m_obj_tree.reset(new KDTree());
	}

	void SPMOEA9::run_() {
		initPop();
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

	void SPMOEA9::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		entry.push_back(m_R1);
		entry.push_back(m_R2);
		entry.push_back(m_R3);
		//Real IGD = m_problem->optimaBase()->invertGenDist(*m_pop);
		entry.push_back(m_IGD);
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA9::updateBuffer() {
		m_solution.clear();
		m_solution.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_solution[0].push_back(&m_pop->at(i));
		/*for (size_t i = 0; i < m_pop->getOffspring().size(); ++i)
			m_solution[1].push_back(&m_pop->getOffspring()[i]);*/
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	void SPMOEA9::initPop() {
		m_pop.reset(new SPMOEA9_pop(m_pop_size, m_problem.get()));
		m_pop->setCR(m_cr);
		m_pop->setMR(m_mr);
		m_pop->setEta(m_ceta, m_meta);

		initiVarSpace(m_problem.get());
		m_pop->initialize(m_problem.get(), m_random.get());
		//时候需要先在子空间进行预采样？

		m_pop->evaluate(m_problem.get(), this);

		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < m_pop->size(); ++i)
			objs.emplace_back(&m_pop->at(i).objective());
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < m_pop->size(); ++i) {
			m_pop->at(i).setFitness(rank[i]);
			if (rank[i] == 0)
				m_history_front_sols.emplace_back(std::make_shared<Solution<>>(m_pop->at(i)));
		}

		//auto& v = *m_param;
		initiObjSpace(m_problem.get());
	}

	SPMOEA9_pop::SPMOEA9_pop(size_t size_pop, Problem *pro): PopSBX<>(size_pop, pro),\
		m_offspring(2 * size_pop, pro, CAST_CONOP(pro)->numberVariables()),
		m_parents(size_pop, pro, CAST_CONOP(pro)->numberVariables())
		{
	}
	
	void SPMOEA9_pop::initialize(Problem *pro, Random *rnd) {
		/*************************/
		/*       种群初始化      */
		/*************************/
		Population<Solution<>>::initialize(pro, rnd);
	}

	void SPMOEA9::initiObjSpace(Problem *pro) {
		for (size_t i = 0; i < m_pop->size(); ++i) {
			m_pop->getOffspring()[i] = m_pop->at(i);
			m_pop->getParent()[i] = m_pop->at(i);
			m_pop->getOffspring()[i + m_pop->size()] = m_pop->at(i);
		}
		/****************************/
		/*      计算目标值范围      */
		/****************************/
		//updateObjRange(pro);
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		m_pop_range.resize(obj_num);
		for (int i = 0; i < obj_num; ++i) {
			m_pop_range[i].first = 1.0e14;
			for (int j = 0; j < m_pop->size(); ++j) {
				if (m_pop->at(j).objective()[i] < m_pop_range[i].first) {
					m_pop_range[i].first = m_pop->at(j).objective()[i];
				}
			}
		}
		for (int i = 0; i < obj_num; ++i) {
			m_pop_range[i].second = -1 * 1.0e14;
			for (int j = 0; j < m_pop->size(); ++j) {
				if (m_pop_range[i].second < m_pop->at(j).objective()[i]) {
					m_pop_range[i].second = m_pop->at(j).objective()[i];
				}
			}
		}
		/****************************/
		/*      更新子目标最优      */
		/****************************/
		for (size_t i = 0; i < obj_num; ++i) {
			m_subobj_opt_sol.emplace_back(m_pop->at(i));
		}
		updateSubObjOpt(m_pop->getParent());
		/****************************/
		/*    目标空间划分初始化    */
		/****************************/
		for (size_t i = 0; i < m_num_region_obj; ++i) {
			m_obj_region_info.emplace_back(new ObjRegionInfo);
		}
		updateObjRange(m_pop->getOffspring(),pro);
		m_front_pop_range = m_pop_range;
		std::vector<std::pair<Real, Real>> m_obj_boundary;
		if (m_normalize) {
			for (size_t i = 0; i < obj_num; ++i) {
				m_obj_boundary.push_back(std::make_pair<>(0., 1.));
			}
		}
		else {
			m_obj_boundary = getFrontPopRange();
		}
		m_mo_hlc->getObjspaceTree().setInitBox(m_obj_boundary);
		m_mo_hlc->getObjspaceTree().inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		m_mo_hlc->getObjspaceTree().buildIndex();
		/*m_obj_tree->setInitBox(m_obj_boundary);
		m_obj_tree->inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		m_obj_tree->buildIndex();*/

		updateObjSpaceInfo(m_pop->getParent(), pro,false);
		std::vector<size_t> space_idx;
		for (size_t i = 0; i < numObjSubspace(); ++i) {
			auto& spaceinfo = getObjRegionInfo(i);
			if (spaceinfo.m_obj_rank == 0) {
				space_idx.push_back(i);
			}
		}
		m_front_space_inx.push_back(space_idx);
	}

	void SPMOEA9::initiVarSpace(Problem *pro) {
		/*************************************/
		/*      MO_HLC搜索空间划分初始化     */
		/*************************************/
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			m_var_boundary.emplace_back(CAST_CONOP(pro)->range(i));
		}
		m_mo_hlc->initialVarSpace(m_var_boundary, m_num_region_var);
		m_pop_var_range = m_var_boundary;
		m_num_e_e.resize(2);
	}

	int SPMOEA9::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		//可视化，策略选项
		auto& v = *m_param;
		m_add_neighbor = v.get<bool>(("cluster add neighbors"));
		m_sample_in_basin = v.get<bool>(("sample in basin"));
		m_subspace_select = v.get<bool>(("select by subspace"));
		m_evolve_by_potential = v.get<bool>(("evolve by potential"));
		/***************************************/
		/*            更新空间信息             */
		/***************************************/
		//使用种群信息更新搜索子空间
		//updateObjSpaceInfo(m_pop->getParent(),pro);
		updateSolSpaceInfo(m_pop->getOffspring(), pro, true);
		updateSolSpaceInfo(m_pop->getParent(),pro,false);
		//对搜索空间的前沿子空间聚类
		auto subspace=spaceAttach(m_pop->getParent(), m_mo_hlc->subspaceTree(), false);
		auto exploit_cluster=clusterExploitSpace(m_mo_hlc->subspaceTree(),subspace[0]);
		//计算剩余空间的探索潜力值，并根据探索潜力值排序
		//getMO_HLC().updateSubspaceInfo();
		//updateVarSpacePotential();
		//对剩余子空间按照探索潜力值聚类
		auto explore_cluster = clusterExploreSpace(getMO_HLC().subspaceTree());
		//更新basin信息，为每个开发潜力区域和探索潜力区域分配搜索潜力和搜索资源
		//m_mo_hlc->initial_basin(exploit_cluster.size()+explore_cluster.size());
		m_mo_hlc->updateBasinInfo(explore_cluster,exploit_cluster,m_mean_potential);
		//根据吸引域搜索潜力为每个吸引域分配搜索资源
		auto potential_inds = predictObjSpacePotential(pro);
		auto resource=assignBasinResource(potential_inds.second);
		
		//E&E ratios
		size_t explore_num = 0;
		std::vector<Real> ee_num;
		for (size_t i = 0; i < explore_cluster.size(); ++i) {
			explore_num += resource[i];
		}
		m_num_e_e[0].push_back((Real)explore_num / m_pop->size());
		m_num_e_e[1].push_back(1-(Real)explore_num / m_pop->size());

		/************************************/
	    /*            性能指标计算          */
		/************************************/
		//recordMetrics(pro);

		/************************************/
		/*            种群演化              */
		/************************************/
		if (this->m_pop->size() % 2 != 0)
			throw MyExcept("population size should be even @NSGAII_SBXRealMu::evolve()");
		int tag = EvaluationTag::kNormalEval;
		//bool b = ifApproxiConverge(rnd);

		/************************************/
		/*            产生子代              */
		/************************************/
		generateOffspring(resource, potential_inds,pro,rnd);

		/***********************************/
		/*            评价子代             */             
		/***********************************/
		//int tag;
		for (size_t i = 0; i < m_pop->size(); i++) {
			tag = m_pop->getOffspring()[i].evaluate(pro,alg);
			auto temp_obj = m_pop->getOffspring()[i].objective();
			for (size_t j = 0; j < temp_obj.size(); ++j) {
				if (temp_obj[j] < 0) {
					size_t a = 1;
				}
			}
			if (tag != EvaluationTag::kNormalEval)
				break;
		}
		/***************************************************/
		/*            更新算法的history和archive           */
		/***************************************************/
		m_pop->NDSort(m_pop->getOffspring());
		
		updateHistorySol(pro);
		Real domain_changed = 0.;
		if (updateHisFrontObjRange()) {
			domain_changed = 1.;
		}

		//updateBuffer();
		/************************************/
		/*             环境选择             */
		/************************************/
		//先更新空间信息
		m_update_tree_flag.push_back(domain_changed);
		updateObjSpaceInfo(m_pop->getOffspring(),pro,domain_changed);//更新子空间信息
		//updateSolSpaceInfo(m_pop->getOffspring(), pro,true);//更新决策空间信息
		updateArchive(m_archive_num,pro);
		//再淘汰选择，采用基于子空间跨度的多参考点选择机制
		envirSelection(pro, rnd);
		//最后更新当前parent和offspring
		for (size_t i = 0; i < m_pop->size(); ++i) {
			m_pop->getOffspring()[i + m_pop->size()] = m_pop->at(i);
			m_pop->getParent()[i] = m_pop->at(i);
		}
		updateObjSpaceInfo(m_pop->getParent(), pro,false);

		//前沿个体索引
		std::vector<std::map<size_t, std::vector<size_t>>> ind=spaceAttach(m_pop->getParent(), m_mo_hlc->getObjspaceTree(), true);
		std::vector<size_t> space_idx;
		for (auto& idx : ind[0]) {
			space_idx.push_back(idx.first);
		}
		m_front_space_inx.push_back(space_idx);

		//记录历史前沿
		//recordHisFront(pro);

		m_pop->iteration()++;
		return tag;
	}

	void SPMOEA9::generateOffspring(const std::vector<size_t>& resource , const std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>& potential_inds, Problem *pro, Random *rnd) {
		/************************************/
		/*            产生新解              */
		/************************************/
		if (!m_evolve_by_potential) {
			//SBX产生新解
			for (size_t i = 0; i < m_pop->size(); i += 2) {
				std::vector<size_t> p(2);
				do {
					p[0] = m_pop->tournamentSelection(pro, rnd);
					p[1] = m_pop->tournamentSelection(pro, rnd);
				} while (p[1] == p[0]);

				//do {//不重复采样
				//	m_pop->crossover(p[0], p[1], m_pop->getOffspring()[i], m_pop->getOffspring()[i + 1], pro, rnd);
				//	m_pop->mutate(m_pop->getOffspring()[i], pro, rnd);
				//	m_pop->mutate(m_pop->getOffspring()[i + 1], pro, rnd);
				//	//子代越界处理
				//	for (size_t j = 0; j < m_pop_var_range.size(); ++j) {
				//		auto range = m_pop_var_range[j];
				//		auto var1 = m_pop->getOffspring()[i].variable().vect()[j];
				//		auto var2 = m_pop->getOffspring()[i + 1].variable().vect()[j];
				//		while (var1 < range.first) {
				//			m_pop->getOffspring()[i].variable().vect()[j] = 2 * range.first - var1;
				//		}
				//		while (var1 > range.second) {
				//			m_pop->getOffspring()[i].variable().vect()[j] = 2 * range.second - var1;
				//		}
				//		while (var2 < range.first) {
				//			m_pop->getOffspring()[i + 1].variable().vect()[j] = 2 * range.first - var1;
				//		}
				//		while (var2 > range.second) {
				//			m_pop->getOffspring()[i + 1].variable().vect()[j] = 2 * range.second - var1;
				//		}
				//	}
				//} while (m_pop->getOffspring()[i].variable().vect() == m_pop->at(p[0]).variable().vect() || \
				//	m_pop->getOffspring()[i].variable().vect() == m_pop->at(p[1]).variable().vect() || \
				//	m_pop->getOffspring()[i + 1].variable().vect() == m_pop->at(p[0]).variable().vect() || \
				//	m_pop->getOffspring()[i + 1].variable().vect() == m_pop->at(p[1]).variable().vect());
				//

				//子代越界处理
				m_pop->crossover(p[0], p[1], m_pop->getOffspring()[i], m_pop->getOffspring()[i + 1], pro, rnd);
				m_pop->mutate(m_pop->getOffspring()[i], pro, rnd);
				m_pop->mutate(m_pop->getOffspring()[i + 1], pro, rnd);
				/*for (size_t j = 0; j < m_pop_var_range.size(); ++j) {
					auto range = m_pop_var_range[j];
					auto var1 = m_pop->getOffspring()[i].variable().vect()[j];
					auto var2 = m_pop->getOffspring()[i + 1].variable().vect()[j];
					while (var1 < range.first) {
						m_pop->getOffspring()[i].variable().vect()[j] = 2 * range.first - var1;
						var1 = m_pop->getOffspring()[i].variable().vect()[j];
					}
					while (var1 > range.second) {
						m_pop->getOffspring()[i].variable().vect()[j] = 2 * range.second - var1;
						var1 = m_pop->getOffspring()[i].variable().vect()[j];
					}
					while (var2 < range.first) {
						m_pop->getOffspring()[i + 1].variable().vect()[j] = 2 * range.first - var1;
						var2=m_pop->getOffspring()[i + 1].variable().vect()[j];
					}
					while (var2 > range.second) {
						m_pop->getOffspring()[i + 1].variable().vect()[j] = 2 * range.second - var1;
						var2 = m_pop->getOffspring()[i + 1].variable().vect()[j];
					}
				}*/
				repairSol(m_pop->getOffspring()[i].variable().vect(),pro,rnd);
				repairSol(m_pop->getOffspring()[i+1].variable().vect(), pro,rnd);
				m_pop->getOffspring()[i].setCounter(0);
				m_pop->getOffspring()[i + 1].setCounter(0);
			}
		}
		else {
			//再根据子空间信息确定采样，按照估计的潜力产生子代（收敛潜力和多样性潜力）
			//个体产生方式是否可以不用GA，按照子空间的潜力值大小分配搜索个体的数目
			//在每个吸引域内产生指定数量的个体：1、基于吸引域内个体生成；2、基于子空间内个体生成
			size_t count_offspring = 0;
			if (m_sample_in_basin) {//在basin内采样
				for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
					//先得到吸引域分配的解个数
					size_t select_num = resource[i];
					auto spaces = m_mo_hlc->getBasinInfo(i).m_subspace_set;
					//再根据根据吸引域的类型和当前解的个数产生新的解
					std::string basin_type = m_mo_hlc->getBasinInfo(i).flag;
					if (basin_type == "explore") {
						//在潜力最大的几个区域随机采样
						std::vector<size_t> temp_space = spaces;
						for (size_t j = 0; j < temp_space.size(); ++j) {
							for (size_t k = j + 1; k < temp_space.size(); ++k) {
								if (m_mo_hlc->getSubspaceInfo(temp_space[j]).m_potential < m_mo_hlc->getSubspaceInfo(temp_space[k]).m_potential) {
									auto temp = temp_space[j];
									temp_space[j] = temp_space[k];
									temp_space[k] = temp;
								}
							}
						}
						size_t count = 0;
						while (count < select_num) {
							for (size_t j = 0; j < temp_space.size(); ++j) {
								if (count < select_num) {
									//在子空间采用随机采样的方式生成一个解
									auto boundary_x = m_mo_hlc->subspaceTree().getBox(temp_space[j]);
									for (size_t k = 0; k < boundary_x.size(); ++k) {
										m_pop->getOffspring()[count_offspring].variable()[k] = boundary_x[k].first + (boundary_x[k].second - boundary_x[k].first) * rnd->uniform.next();
									}
									m_pop->getOffspring()[count_offspring].setId(count_offspring);
									m_pop->getOffspring()[count_offspring].setCounter(0);
									count_offspring++;
									count++;
								}
								else {
									break;
								}
							}
						}
					}
					else {
						//使用吸引域内部当前个体产生解
						size_t ind_num = m_mo_hlc->getBasinInfo(i).m_current_indi.size();
						/*****************采用SBX产生*******************/
						PopSBX<> temp_pop(ind_num, pro);
						//初始化种群
						size_t count = 0;
						for (size_t k = 0; k < spaces.size(); ++k) {
							size_t num = m_mo_hlc->getSubspaceInfo(spaces[k]).m_curr_sols.size();
							if (num != 0) {
								for (size_t p = 0; p < num; ++p) {
									temp_pop.at(count) = *m_mo_hlc->getSubspaceInfo(spaces[k]).m_curr_sols[p];
									count++;
								}
							}
						}
						if (ind_num > 1) {
							for (size_t j = 0; j < select_num; ++j) {
								std::vector<size_t> p(ind_num, 0);
								std::vector<Real> var1;
								std::vector<Real> var2;
								Real dist = 0.;
								do {
									if (ind_num > 3) {
										p[0] = temp_pop.tournamentSelection(pro, rnd);
										p[1] = temp_pop.tournamentSelection(pro, rnd);
									}
									else {
										for (size_t j = 0; j < ind_num; ++j) {
											p[j] = j;
										}
										rnd->uniform.shuffle(p.begin(), p.end());
									}
									if (p[0] != p[1]) {
										var1 = temp_pop[p[0]].variable().vect();
										var2 = temp_pop[p[1]].variable().vect();
										dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
										if (dist < 1. * 10e-6) {
											size_t h_idx = std::floor(m_mo_hlc->getBasinInfo(i).m_history_inds.size() * rnd->uniform.next());
											temp_pop.append(*m_mo_hlc->getBasinInfo(i).m_history_inds[h_idx]);
											p.resize(ind_num + 1);
											ind_num++;
										}
									}
								} while (p[1] == p[0] || dist <= 1. * 10e-6);
								//判断交叉之后子代时候与父代相同，相同则继续
								do {
									temp_pop.crossover(p[0], p[1], m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
									temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
									//子代越界处理
									//for (size_t j = 0; j < m_pop_var_range.size(); ++j) {
									//	auto range = m_pop_var_range[j];
									//	auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
									//	while (var1 < range.first) {
									//		m_pop->getOffspring()[count_offspring].variable().vect()[j] = 2 * range.first - var1;
									//		var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
									//	}
									//	while (var1 > range.second) {
									//		m_pop->getOffspring()[count_offspring].variable().vect()[j] = 2 * range.second - var1;
									//		var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
									//	}
									//	//判断是否越界
									//	auto tmp = m_pop->getOffspring()[count_offspring].variable().vect()[j];
									//	int isnan(tmp);
									//	if (!isnormal(tmp)) {
									//		size_t a = 1;
									//	}
									//	if (tmp<range.first || tmp>range.second) {
									//		size_t a = 1;
									//	}
									//}
									repairSol(m_pop->getOffspring()[count_offspring].variable().vect(),pro,rnd);
								} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[0]].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[1]].variable().vect());

								m_pop->getOffspring()[count_offspring].setId(count_offspring);
								m_pop->getOffspring()[count_offspring].setCounter(0);
								count_offspring++;
							}
						}
						else {
							//采用局部变异的方式，先找出个体所在子空间
							size_t idx = 0;
							for (size_t j = 0; j < spaces.size(); ++j) {
								if (!m_mo_hlc->getSubspaceInfo(spaces[j]).m_curr_sols.empty()) {
									idx = j;
								}
							}
							auto temp_x = m_mo_hlc->getSubspaceInfo(spaces[idx]).m_curr_sols.back()->variable().vect();
							auto boundary_x = m_mo_hlc->subspaceTree().getBox(spaces[idx]);
							for (size_t j = 0; j < select_num; ++j) {
								auto temp_off = localMutationSol(temp_x, boundary_x,pro, rnd);
								m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
								m_pop->getOffspring()[count_offspring].setId(count_offspring);
								m_pop->getOffspring()[count_offspring].setCounter(0);
								count_offspring++;
							}
						}
					}
				}
			}
			else {//在子空间采样
				for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
					//在吸引域中，根据子空间的搜索潜力决定在哪个子空间采样
					std::string basin_type = m_mo_hlc->getBasinInfo(i).flag;
					size_t select_num = resource[i];
					std::vector<size_t> space_idx = m_mo_hlc->getBasinInfo(i).m_subspace_set;
					std::vector<Real> subspace_potential;
					Real sum = 0.;
					for (size_t j = 0; j < space_idx.size(); ++j) {
						subspace_potential.push_back(m_mo_hlc->getSubspaceInfo(space_idx[j]).m_potential);
						sum += subspace_potential.back();
					}
					for (auto& i : subspace_potential) {
						i = i / sum;
					}
					auto select_probability = subspace_potential;
					for (size_t i = 0; i < subspace_potential.size(); ++i) {
						Real sum = 0.;
						for (size_t j = 0; j <= i; ++j) {
							sum += subspace_potential[j];
						}
						select_probability[i] = sum;
					}
					for (size_t i = 0; i < select_num; ++i) {
						Real rand_num = rnd->uniform.next();
						size_t idx;
						for (size_t j = 0; j < select_probability.size(); ++j) {
							if (rand_num <= select_probability[j]) {
								idx = j;
								break;
							}
						}
						//选择采样的子空间为space_idx[idx],并在里面进行采样
						auto boundary_x = m_mo_hlc->subspaceTree().getBox(space_idx[idx]);
						if (basin_type == "explore" || m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.empty()) {//随机采样的方式
							for (size_t k = 0; k < boundary_x.size(); ++k) {
								m_pop->getOffspring()[count_offspring].variable()[k] = boundary_x[k].first + (boundary_x[k].second - boundary_x[k].first) * rnd->uniform.next();
							}
						}
						else {//基于历史解生成，有问题，生成解
							if (m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.size() > 1) {
								size_t pop_num = m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.size();
								/*****************采用SBX产生*******************/
								PopSBX<> temp_pop(pop_num, pro);
								//初始化种群
								for (size_t j = 0; j < pop_num; ++j) {
									temp_pop.at(j) = *m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols[j];
								}
								std::vector<size_t> p(pop_num, 0);
								std::vector<Real> var1;
								std::vector<Real> var2;
								Real dist = 0.;
								do {
									if (pop_num > 3) {
										p[0] = temp_pop.tournamentSelection(pro, rnd);
										p[1] = temp_pop.tournamentSelection(pro, rnd);
									}
									else {
										for (size_t j = 0; j < pop_num; ++j) {
											p[j] = j;
										}
										rnd->uniform.shuffle(p.begin(), p.end());
									}
									if (p[0] != p[1]) {
										var1 = temp_pop[p[0]].variable().vect();
										var2 = temp_pop[p[1]].variable().vect();
										dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
										if (dist < 1. * 10e-6) {
											size_t h_idx = std::floor(m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds.size() * rnd->uniform.next());
											temp_pop.append(*m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds[h_idx]);
											p.resize(pop_num + 1);
											pop_num++;
										}
									}
								} while (p[1] == p[0] || dist <= 1. * 10e-6);
								do {
									temp_pop.crossover(p[0], p[1], m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
									temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
									//子代越界处理
									/*for (size_t k = 0; k < m_pop_var_range.size(); ++k) {
										auto range = m_pop_var_range[k];
										auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
										while (var1 < range.first) {
											m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.first - var1;
											var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
										}
										while (var1 > range.second) {
											m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.second - var1;
											auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
										}
									}*/
									repairSol(m_pop->getOffspring()[count_offspring].variable().vect(),pro,rnd);
								} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[0]].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[1]].variable().vect());

								/*****************采用DE产生*******************/
								//PopDE<> temp_pop(pop_num,pro);
								////初始化种群
								////交叉变异

							}
							else {//当只有一个解时，在解的邻域产生,每一维的3 sigma
								auto temp_x = m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.back()->variable().vect();
								auto boundary_x = m_mo_hlc->subspaceTree().getBox(space_idx[idx]);
								auto temp_off = localMutationSol(temp_x, boundary_x,pro, rnd);
								m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
							}
						}
						m_pop->getOffspring()[count_offspring].setId(count_offspring);
						m_pop->getOffspring()[count_offspring].setCounter(0);
						count_offspring++;
					}
				}
			}
			//在稀疏点处采样
			size_t num_basin = m_mo_hlc->numBasin();
			size_t num_obj = pro->numberObjectives();
			std::vector<std::vector<size_t>> poten_inds = potential_inds.second;
			for (size_t i =num_basin ; i < resource.size();++i) {
				if (i >= num_basin&&i<num_basin+num_obj) {
					//依次得到每个点所在的子空间，然后向延伸方向采样
					auto pop_idx = poten_inds[i - num_basin][0];
					auto var_position = m_pop->at(pop_idx).variable().vect();
					auto obj_position = m_pop->at(pop_idx).objective();
					auto var_index = m_mo_hlc->subspaceTree().getRegionIdx(var_position);
					auto boundary_x = m_mo_hlc->subspaceTree().getBox(var_index);
					auto obj_index =m_mo_hlc->getObjspaceTree().getRegionIdx(obj_position);
					auto ind_index = getObjRegionInfo(obj_index).ind_idx;
					if (ind_index.size() > 1) {
						std::vector<Real> obj1= m_pop->at(pop_idx).objective();;
						std::vector<Real> obj2= m_pop->at(ind_index[0]).objective();;
						Real dist = 0.;
						size_t count = 0;
						std::vector<Real> ind_var;
						do {
							rnd->uniform.shuffle(ind_index.begin(), ind_index.end());
							obj2 = m_pop->at(ind_index[0]).objective();
							dist = euclideanDistance(obj1.begin(), obj1.end(), obj2.begin());
							if (dist < 1. * 10e-6) {
								count++;
							}
							if (count > 20) {
								break;
							}
						} while (ind_index[0] == pop_idx|| dist < 1. * 10e-7);
						ind_var = m_pop->at(ind_index[0]).variable().vect();
						if (count > 20) {
							size_t his_sol_num = m_mo_hlc->getSubspaceInfo(var_index).m_history_inds.size();
							std::vector<Real> var1;
							do {
								size_t select_idx = std::floor(his_sol_num * rnd->uniform.next());
								var1 = (*m_mo_hlc->getSubspaceInfo(var_index).m_history_inds[select_idx]).variable().vect();
								dist = euclideanDistance(var1.begin(), var1.end(), var_position.begin());
							} while (dist < 1. * 10e-4);
							ind_var = var1;
						}
						auto temp_off = vectorMutationSol(var_position, ind_var,pro, rnd);
						m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
						m_pop->getOffspring()[count_offspring].setId(count_offspring);
						m_pop->getOffspring()[count_offspring].setCounter(0);
						count_offspring++;
					}
					//if (m_mo_hlc->getSubspaceInfo(var_index).m_curr_sols.size() > 1) {
					//	//方向采样
					//	auto ind_idx = m_mo_hlc->getSubspaceInfo(var_index).m_curr_sols_idx;
					//	for (size_t k = 0; k < resource[i]; ++k) {
					//		//任意找一个非边界点,但不和边界点重合
					//		std::vector<Real> var1;
					//		std::vector<Real> var2;
					//		do {
					//			rnd->uniform.shuffle(ind_idx.begin(), ind_idx.end());
					//			var1 = temp_pop[p[0]].variable().vect();
					//			var2 = temp_pop[p[1]].variable().vect();
					//			dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
					//			if (dist < 1. * 10e-6) {
					//				size_t h_idx = std::floor(m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds.size() * rnd->uniform.next());
					//				temp_pop.append(*m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds[h_idx]);
					//				p.resize(pop_num + 1);
					//				pop_num++;
					//			}
					//		} while (ind_idx[0] == pop_idx);
					//		//向量计算，步长设置
					//		auto select_ind = m_pop->at(ind_idx[0]).variable().vect();
					//		auto temp_off = vectorMutationSol(var_position,select_ind,rnd);
					//		m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
					//		m_pop->getOffspring()[count_offspring].setId(count_offspring);
					//		m_pop->getOffspring()[count_offspring].setCounter(0);
					//		count_offspring++;
					//	}
					//}
					else {
						//局部变异
						for (size_t k = 0; k < resource[i]; ++k) {
							auto temp_off = localMutationSol(var_position, boundary_x,pro, rnd);
							m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
							m_pop->getOffspring()[count_offspring].setId(count_offspring);
							m_pop->getOffspring()[count_offspring].setCounter(0);
							count_offspring++;
						}
					}
				}
				else {
					//先得到两个个体及其所在的子空间
					//auto ind1 = m_pop->at((*poten_inds[i-num_basin+num_obj])[0]).variable().vect();
					//auto ind2 = m_pop->at((*poten_inds[i - num_basin + num_obj])[1]).variable().vect();
					/*****************采用SBX产生*******************/
					PopSBX<> temp_pop(2, pro);
					temp_pop.setEta(2,2);
					//初始化种群
					for (size_t j = 0; j < 2; ++j) {
						temp_pop.at(j) = m_pop->at(poten_inds[i-num_basin][j]);
					}
					for (size_t j = 0; j < resource[i]; ++j) {
						do {
							temp_pop.crossover(0, 1, m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
							temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
							//子代越界处理
							for (size_t k = 0; k < m_pop_var_range.size(); ++k) {
								auto range = m_pop_var_range[k];
								auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
								while (var1 < range.first) {
									m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.first - var1;
									var1=m_pop->getOffspring()[count_offspring].variable().vect()[k];
								}
								while (var1 > range.second) {
									m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.second - var1;
									var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
								}
								////越界判断
								//auto tmp = m_pop->getOffspring()[count_offspring].variable().vect()[k];
								//if (!isnormal(tmp)) {
								//	size_t a = 1;
								//}
								//if (tmp<m_pop_var_range[k].first || tmp>m_pop_var_range[k].second) {
								//	size_t a = 1;
								//}
							}
							repairSol(m_pop->getOffspring()[count_offspring].variable().vect(),pro,rnd);
						} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[0].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[1].variable().vect());
						m_pop->getOffspring()[count_offspring].setId(count_offspring);
						m_pop->getOffspring()[count_offspring].setCounter(0);
						count_offspring++;
					}
				}
			}
		}
	}

	std::vector<Real> SPMOEA9::localMutationSol(const std::vector<Real>& sol, const std::vector<std::pair<Real,Real>>& boundary,Problem *pro, Random *rnd) {
		std::vector<std::pair<Real, Real>> dim_dist;
		for (size_t j = 0; j < sol.size(); ++j) {
			auto temp1 = (sol[j] - boundary[j].first) / 3;
			auto temp2 = (boundary[j].second - sol[j]) / 3;
			dim_dist.push_back(std::make_pair<>(temp1, temp2));
		}
		std::vector<Real> inds;
		for (size_t k = 0; k < boundary.size(); ++k) {
			Real tmp = boundary[k].first;
			while ((tmp - boundary[k].first) * (tmp - boundary[k].second) >= 0) {
				Real rand_num = rnd->uniform.next();
				if (rand_num > 0.5) {
					tmp = rnd->normal.nextNonStd(sol[k], dim_dist[k].second);
					if (tmp < sol[k]) {
						tmp = 2 * sol[k] - tmp;
					}
				}
				else {
					tmp = rnd->normal.nextNonStd(sol[k], dim_dist[k].first);
					if (tmp > sol[k]) {
						tmp = 2 * sol[k] - tmp;
					}
				}
			}
			////判断是否越界
			//if (!isnormal(tmp)) {
			//	size_t a = 1;
			//}
			//if (tmp<m_pop_var_range[k].first || tmp>m_pop_var_range[k].second) {
			//	size_t a = 1;
			//}
			inds.push_back(tmp);
		}
		repairSol(inds,pro,rnd);
		return inds;
	}

	std::vector<Real> SPMOEA9::vectorMutationSol(const std::vector<Real>& sol1, const std::vector<Real>& sol2, Problem *pro, Random *rnd) {
		//先求穿出点坐标
		std::vector<Real> boundary_points;
		for (size_t i = 0; i < sol1.size(); ++i) {
			boundary_points.clear();
			Real temp_ratio1 = (m_pop_var_range[i].first - sol1[i]) / (sol1[i]- sol2[i]);
			if (temp_ratio1 >= 0) {
				for (size_t j = 0; j < sol1.size(); ++j) {
					if (j == i) {
						boundary_points.push_back(m_pop_var_range[i].first);
					}
					else {
						auto temp = sol1[j] + temp_ratio1 * (sol1[j] - sol2[j]);
						if (temp<m_pop_var_range[j].first || temp>m_pop_var_range[j].second) {
							break;
						}
						else {
							boundary_points.push_back(temp);
						}
					}
				}
			}
			if (boundary_points.size() == sol1.size()) {
				break;
			}
			boundary_points.clear();
			Real temp_ratio2 = (m_pop_var_range[i].second - sol1[i]) / (sol1[i] - sol2[i]);
			if (temp_ratio2 >= 0) {
				for (size_t j = 0; j < sol1.size(); ++j) {
					if (j == i) {
						boundary_points.push_back(m_pop_var_range[i].second);
					}
					else {
						auto temp = sol1[j] + temp_ratio2 * (sol1[j] - sol2[j]);
						if (temp<m_pop_var_range[j].first || temp>m_pop_var_range[j].second) {
							break;
						}
						else {
							boundary_points.push_back(temp);
						}
					}
				}
			}
			if (boundary_points.size() == sol1.size()) {
				break;
			}
		}
		//再由穿出点与边界点随机合成
		Real temp_rand = rnd->uniform.next();
		std::vector<Real> temp_off;
		for (size_t i = 0; i < sol1.size(); ++i) {
			temp_off.push_back(sol1[i]+temp_rand*(boundary_points[i]-sol1[i]));
		}
		repairSol(temp_off,pro,rnd);
		return temp_off;
	}

	void SPMOEA9::repairSol(std::vector<Real>& sol, Problem *pro, Random *rnd) {
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			m_var_boundary.emplace_back(CAST_CONOP(pro)->range(i));
		}
		for (size_t i = 0; i < sol.size(); ++i) {
			if (sol[i] < m_var_boundary[i].first || (!isnormal(sol[i]))) {
				sol[i] = m_var_boundary[i].first + 0.5 * rnd->uniform.next() * (m_var_boundary[i].second - m_var_boundary[i].first);
			}
			if (sol[i] > m_var_boundary[i].second || (!isnormal(sol[i]))) {
				sol[i] = m_var_boundary[i].second - 0.5 * rnd->uniform.next() * (m_var_boundary[i].second - m_var_boundary[i].first);
			}

			/*if (sol[i] < m_pop_var_range[i].first || (!isnormal(sol[i]))) {
				sol[i] = m_pop_var_range[i].first + 0.5 * rnd->uniform.next() * (m_pop_var_range[i].second - m_pop_var_range[i].first);
			}
			if (sol[i] > m_pop_var_range[i].second || (!isnormal(sol[i]))) {
				sol[i] = m_pop_var_range[i].second - 0.5 * rnd->uniform.next() * (m_pop_var_range[i].second - m_pop_var_range[i].first);
			}*/
		}
	}

	void SPMOEA9::updateObjSpaceInfo(const Population<Solution<>>& pop,Problem *pro,bool changed) {
		/*更新目标值范围*/
		updateSubObjOpt(pop);
		bool b=updateObjRange(pop,pro);
		/*更新目标空间划分*/
		if (changed&&!m_normalize) {
			updateObjTree();
		}
		/**********************************/
		/*      目标子空间加入个体信息    */
		/**********************************/
		/*目标空间清空,搜索空间选择性清空*/
		m_obj_region_info.clear();
		for (size_t i = 0; i < m_num_region_obj; ++i) {
			m_obj_region_info.emplace_back(new ObjRegionInfo);
		}
		/*先更新目标子空间*/
		for (size_t i = 0; i < pop.size(); i++) {
			auto temp_obj = pop[i].objective();
			if (m_normalize) {
				for (size_t j = 0; j < temp_obj.size(); ++j) {
					temp_obj[j] = (temp_obj[j]- m_front_pop_range[j].first) / (m_front_pop_range[j].second- m_front_pop_range[j].first);
				}
			}
			size_t idx_obj_region = m_mo_hlc->getObjspaceTree().getRegionIdx(temp_obj);
			if (m_obj_region_info[idx_obj_region]->obj_optima.empty()) {//obj
				for (size_t j = 0; j < temp_obj.size(); j++) {
					m_obj_region_info[idx_obj_region]->obj_optima.emplace_back(pop[i].objective()[j]);
				}
			}
			else {
				for (size_t j = 0; j < temp_obj.size(); j++) {
					Real temp = m_obj_region_info[idx_obj_region]->obj_optima[j];
					m_obj_region_info[idx_obj_region]->obj_optima[j]=std::min(temp, pop[i].objective()[j]);
				}
			}
			m_obj_region_info[idx_obj_region]->m_curr_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));//solution
		    //个体的索引
			m_obj_region_info[idx_obj_region]->ind_idx.push_back(i);																								 
			if (pop[i].fitness() < m_obj_region_info[idx_obj_region]->m_obj_rank) { //目标子空间排序
				m_obj_region_info[idx_obj_region]->m_obj_rank = pop[i].fitness();
			}
		}
	}

	void SPMOEA9::updateSubObjOpt(const Population<Solution<>>& pop) {
		for (size_t i = 0; i < m_subobj_opt_sol.size(); ++i) {
			size_t inx = 0;
			Real temp_obj = m_subobj_opt_sol[i].objective()[i];
			for (size_t j = 0; j < pop.size(); ++j) {//第i个解存放第i个目标的最小值
				if (temp_obj > pop[j].objective()[i]) {
					temp_obj = pop[j].objective()[i];
					inx = j;
				}
			}
			if (m_subobj_opt_sol[i].objective()[i] > pop[inx].objective()[i]) {
				m_subobj_opt_sol[i] = pop[inx];
			}
		}
	}

	bool SPMOEA9::updateObjRange(const Population<Solution<>>& pop, Problem *pro) {
		std::vector<std::pair<Real, Real>> pop_range_before = getPopRange();
		int M = CAST_CONOP(pro)->numberObjectives();
		int N = CAST_CONOP(pro)->numberVariables();
		if (getPopRange().empty()) {
			m_pop_range.resize(M);
		}
		//更新当前种群范围
		for (int i = 0; i < M; ++i) {
			m_pop_range[i].first = 1.0e14;
			m_pop_range[i].second = -1*1.0e14;
			for (int j = 0; j < pop.size(); ++j) {
				if (m_pop_range[i].first > pop[j].objective()[i]) {
					m_pop_range[i].first = pop[j].objective()[i];
				}
				if (m_pop_range[i].second < pop[j].objective()[i]) {
					m_pop_range[i].second = pop[j].objective()[i];
				}
			}
		}
		////采用子目标历史最优和当前种群范围更新
		//for (size_t i = 0; i < M; ++i) {
		//	m_pop_range[i].first = m_subobj_opt_sol[i].objective()[i];
		//}
		//for (int i = 0; i < M; ++i) {
		//	m_pop_range[i].second = -1 * 1.0e14;
		//	for (int j = 0; j < pop.size(); ++j) {
		//		if (m_pop_range[i].second < pop[j].objective()[i]) {
		//			m_pop_range[i].second = pop[j].objective()[i];
		//		}
		//	}
		//}
		////archive张成的空间范围
		//for (int i = 0; i < M; ++i) {
		//	m_pop_range[i].first = 1.0e14;
		//	for (int j = 0; j < m_history_front_sols.size(); ++j) {
		//		if (m_history_front_sols[j]->objective()[i] < m_pop_range[i].first) {
		//			m_pop_range[i].first = m_history_front_sols[j]->objective()[i];
		//		}
		//	}
		//}
		//for (int i = 0; i < M; ++i) {
		//	m_pop_range[i].second = -1 * 1.0e14;
		//	for (int j = 0; j < m_history_front_sols.size(); ++j) {
		//		if (m_pop_range[i].second < m_history_front_sols[j]->objective()[i]) {
		//			m_pop_range[i].second = m_history_front_sols[j]->objective()[i];
		//		}
		//	}
		//}
		////使用历史子目标的极值点更新范围
		//std::vector<Solution<>> temp_pop;
		//for (size_t i = 0; i < m_pop->getOffspring().size(); ++i) {
		//	temp_pop.emplace_back(m_pop->getOffspring()[i]);
		//}
		//
		//if (getPopRange().empty()) {
		//	m_pop_range.resize(M);
		//}
		//int size = temp_pop.size();
		////取历史各子目标的最值点，然后计算边界
		//for (int i = 0; i < M; ++i) {
		//	//m_pop_range[i].first = 1.0e14;
		//	for (int j = 0; j < size; ++j) {
		//		if (temp_pop[j].objective()[i] < m_pop_range[i].first) {
		//			m_pop_range[i].first = temp_pop[j].objective()[i];
		//		}
		//	}
		//}
		//for (int i = 0; i < M; ++i) {
		//	m_pop_range[i].second = -1 * 1.0e14;
		//	for (int j = 0; j < size; ++j) {
		//		if (m_pop_range[i].second < temp_pop[j].objective()[i]) {
		//			m_pop_range[i].second = temp_pop[j].objective()[i];
		//		}
		//	}
		//}
		bool changed = false;
		if (empty(pop_range_before)) {
			return true;//范围发生了变化
		}
		else {
			for (size_t i = 0; i < pop_range_before.size(); ++i) {
				if (pop_range_before[i].first != m_pop_range[i].first || pop_range_before[i].second != m_pop_range[i].second)
					return true;
			}
			return false;
		}
	}

	void SPMOEA9::updateObjTree() {
		//auto obj_range = getPopRange();
		auto obj_range = getFrontPopRange();
		m_mo_hlc->getObjspaceTree().setInitBox(obj_range);
		m_mo_hlc->getObjspaceTree().inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		m_mo_hlc->getObjspaceTree().buildIndex();
	}

	void SPMOEA9::updateSolSpaceInfo(const Population<Solution<>>& pop, Problem *pro,bool b) {
		//没有解和其他信息的子空间，能否由邻域信息推断？
		//得到优化模式
		//依次遍历解来更新子空间
		//有些信息需要及时清空
		if (b) {//b=true表示更新历史信息
			size_t num = pop.size() / 2;
			for (size_t i = 0; i < num; ++i) {
				auto temp_var = pop[i].variable().vect();
				auto temp_obj = pop[i].objective();
				size_t idx_var_region = m_mo_hlc->subspaceTree().getRegionIdx(temp_var);
				//加入子目标最优值
				if (m_mo_hlc->getSubspaceInfo(idx_var_region).m_history_inds.empty()) {//var
					for (size_t j = 0; j < temp_obj.size(); j++) {
						m_mo_hlc->getSubspaceInfo(idx_var_region).m_subObj_optima.push_back(temp_obj[j]);
						//m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve.push_back(0.);
					}
				}
				else {
					for (size_t j = 0; j < temp_obj.size(); j++) {
						Real temp = m_mo_hlc->getSubspaceInfo(idx_var_region).m_subObj_optima[j];
						if (temp > temp_obj[j]) {
							m_mo_hlc->getSubspaceInfo(idx_var_region).m_subObj_optima[j] = temp_obj[j];
							//Real temp_improve = m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve[j];
							//m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve[j] = std::min(temp_improve,temp_obj[j]-temp);
						}
					}
				}
				//update historical sols
				m_mo_hlc->getSubspaceInfo(idx_var_region).m_history_inds.emplace_back(std::make_shared<Solution<>>(pop[i]));
				//update sample frequency
				m_mo_hlc->getSubspaceInfo(idx_var_region).m_sub_freq++;
			}
			//加入每代的前沿个体
			Population<Solution<>> temp_pop;
			for (size_t i = 0; i < num; ++i) {
				if (pop[i].fitness() == 0) {
					temp_pop.append(pop[i]);
				}
			}
			m_gen_front_sols.emplace_back(std::make_shared<Population<Solution<>>>(temp_pop));
		}
		else {//更新当前信息
			for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
				m_mo_hlc->getSubspaceInfo(i).m_curr_sols.clear();
				//m_mo_hlc->getSubspaceInfo(i).m_curr_sols_idx.clear();
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = INT16_MAX;
			}
			int max_rank = 0;
			for (size_t i = 0; i < pop.size(); i++) {
				auto temp_var = pop[i].variable().vect();
				auto temp_obj = pop[i].objective();
				size_t idx_var_region = m_mo_hlc->subspaceTree().getRegionIdx(temp_var);
				//加入子空间内的最好排序值
				if (m_mo_hlc->getSubspaceInfo(idx_var_region).m_curr_sols.empty()) {//rank
					m_mo_hlc->getSubspaceInfo(idx_var_region).m_best_rank = pop[i].fitness();
					if (pop[i].fitness() > max_rank) {
						max_rank = pop[i].fitness();
					}
				}
				else {
					m_mo_hlc->getSubspaceInfo(idx_var_region).m_best_rank = std::min(m_mo_hlc->getSubspaceInfo(idx_var_region).m_best_rank, pop[i].fitness());
					if (m_mo_hlc->getSubspaceInfo(idx_var_region).m_best_rank > max_rank) {
						max_rank = m_mo_hlc->getSubspaceInfo(idx_var_region).m_best_rank;
					}
				}
				//加入子空间内的当前解
				m_mo_hlc->getSubspaceInfo(idx_var_region).m_curr_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));//solution
				//m_mo_hlc->getSubspaceInfo(idx_var_region).m_curr_sols_idx.push_back(i);
			    //加入子空间内的前沿解
				if (pop[i].fitness() == 0) {
					m_mo_hlc->getSubspaceInfo(idx_var_region).m_front_sol_sub.push_back(std::make_shared<Solution<>>(pop[i]));
				}
			}

			//为没有当前个体的子空间分配rank值
			for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
				auto& spaceinfo = m_mo_hlc->getSubspaceInfo(i);
				if (spaceinfo.m_best_rank == INT16_MAX) {
					spaceinfo.m_best_rank = max_rank + 1;
				}
			}

			////更新探索子空间中子目标的排序值
			//std::vector<size_t> explore_space = getExploreSpace();
			////updateSubObjRank(explore_space, pro);
			////更新子空间中子目标的排序值
			//std::vector<size_t> exploit_space = getExploitSpace();
			//std::vector<size_t> spaces = explore_space;
			//for (auto& i : exploit_space) {
			//	spaces.push_back(i);
			//}
			//更新子空间中子目标的排序值
			updateSubObjRank(pro);
			//更新子空间潜力,先根据目标空间分布，找出有潜力的子空间(边界就空隙大的区域)，再得到这些潜力子空间的中间区域
			auto obj_var_potential = predictObjSpacePotential(pro);//目标空间边界点和稀疏点
			calSpacePotential(obj_var_potential.first,pro);
			/*for (size_t i = 0; i < explore_space.size(); ++i) {
				m_mo_hlc->getSubspaceInfo(explore_space[i]).m_explore_potential = m_mo_hlc->calExploreSpacePotential(explore_space[i]);
			}
			for (size_t i = 0; i < exploit_space.size(); ++i) {
				m_mo_hlc->getSubspaceInfo(exploit_space[i]).m_exploit_potential = m_mo_hlc->calExploitSpacePotential(exploit_space[i], obj_var_potential);
			}*/
		}
	}

	void SPMOEA9_pop::NDSort(Population<Solution<>>& pop) {
		nondominatedSorting(pop);
	}

	void SPMOEA9::updateArchive(size_t num, Problem *pro) {
		//从历史前沿中选择一定数量的个体作为archive
		//三种选择方式：1、拥挤选择；2、参考向量选择；3、子空间跨度选择
		std::vector<size_t> select_index;
	}

	void SPMOEA9::updateHistorySol(Problem *pro) {
		std::vector<Solution<>> temp_pop;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			if (m_pop->getOffspring()[i].fitness() == 0) {
				temp_pop.emplace_back(m_pop->getOffspring()[i]);
			}
		}
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

	bool SPMOEA9::updateHisFrontObjRange() {
		if (m_history_front_sols.size() > 1) {
			auto temp_range = m_front_pop_range;
			m_front_pop_range.clear();
			bool changed = false;
			for (size_t i = 0; i < temp_range.size(); ++i) {
				std::pair<Real, Real> dim_range;
				Real min_value = 1. * 10e14;
				Real max_value = -1. * 10e14;
				for (size_t j = 0; j < m_history_front_sols.size(); ++j) {
					auto temp_obj = m_history_front_sols[j]->objective();
					if (temp_obj[i] < min_value) {
						min_value = temp_obj[i];
					}
					if (temp_obj[i] > max_value) {
						max_value = temp_obj[i];
					}
				}
				dim_range.first = min_value;
				dim_range.second = max_value;
				m_front_pop_range.emplace_back(dim_range);
			}
			for (size_t i = 0; i < temp_range.size(); ++i) {
				if (temp_range[i].first != m_front_pop_range[i].first || temp_range[i].second != m_front_pop_range[i].second) {
					changed = true;
					break;
				}
			}
			return changed;
		}
		else {
			return false;
		}
	}

	void SPMOEA9::calSpacePotential(std::vector<size_t>& spaces,Problem *pro) {
		//根据子空间的采样反馈系数、最佳排序值、子目标的最好排序值,以及目标空间的特征影响的潜力值
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			auto& spaceinfo = m_mo_hlc->getSubspaceInfo(i);
			auto sample_feedback = spaceinfo.m_coff_feedback;
			
			size_t nd_rank = spaceinfo.m_best_rank;
			auto subobjrank = spaceinfo.m_subObj_optima_rank;//目标值越好，rank=1
			size_t best_rank = INT16_MAX;
			for (size_t i = 0; i < subobjrank.size(); ++i) {
				if (best_rank > subobjrank[i]) {
					best_rank = subobjrank[i];
				}
			}
			int fre = spaceinfo.m_sub_freq;
			Real potential=1./std::pow(1+nd_rank,0.5)+1./std::pow(best_rank,0.1);
			Real num = 30*(pro->numberVariables()-1);
			Real a = 10,b=1,c=10./num,d=7;
			/*if (std::find(spaces.begin(), spaces.end(), i) != spaces.end()) {
				potential += 0.5;
			}*/
			spaceinfo.m_exploit_potential = potential;
			Real explore_potential = a / (b + std::exp(c * fre - d));
			spaceinfo.m_explore_potential = explore_potential;
			potential += explore_potential;
			potential *= sample_feedback;
			spaceinfo.m_potential = potential;
		}
	}

	void SPMOEA9::updateSubObjRank(Problem *pro) {
		//先得到这些子空间目标值的最差值
		std::vector<Real> worst_value;
		//然后为各子空间排序
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		std::vector<size_t> has_counted(obj_num);
		for (size_t i = 0; i < obj_num; ++i) {
			size_t count_sort = 0;
			for (size_t j = 0; j < m_mo_hlc->numSubspace(); ++j) {
				Real temp_obj;
				m_mo_hlc->getSubspaceInfo(j).m_subObj_optima_rank.resize(obj_num);
				if (m_mo_hlc->getSubspaceInfo(j).m_subObj_optima.empty()) {
					continue;
				}
				else {
					temp_obj = m_mo_hlc->getSubspaceInfo(j).m_subObj_optima[i];
					count_sort++;
				}
				size_t count = 1;
				for (size_t k = 0; k < m_mo_hlc->numSubspace(); ++k) {
					if (k!= j) {
						if (m_mo_hlc->getSubspaceInfo(k).m_subObj_optima.empty()) {
							continue;
						}
						else if(temp_obj> m_mo_hlc->getSubspaceInfo(k).m_subObj_optima[i]){
							count++;
						}
					}
				}
				m_mo_hlc->getSubspaceInfo(j).m_subObj_optima_rank[i] = count;
			}
			has_counted[i] = count_sort;
		}
		for (size_t i = 0; i < obj_num; ++i) {
			for (size_t j = 0; j < m_mo_hlc->numSubspace(); ++j) {
				if (m_mo_hlc->getSubspaceInfo(j).m_subObj_optima.empty()) {
					m_mo_hlc->getSubspaceInfo(j).m_subObj_optima_rank[i] = has_counted[i] + 1;
				}
			}
		}
	}

	bool SPMOEA9::ifApproxiConverge(Random *rnd) {
		//通过不同层个体的比例来反映
		//产生一个随机数与m_R1比较
		Real rand_num = rnd->uniform.next();
		if (rand_num > m_R1) {
			return true;
		}
		else {
			return false;
		}
		//size_t count = 0;
		//for (size_t i = 0; i < m_pop->size(); i++) {
		//	if (m_pop->at(i).surviveAge() >= m_converge_age)//连续存活n代
		//		count++;
		//}
		//if (count >= m_pop->size() / 10)//连续存活的个体数达到一定规模
		//	return true;
		//else
		//	return false;
	}

	bool SPMOEA9::popConverged() {
		size_t count = 0;
		for (size_t i = 0; i < m_pop->getOffspring().size(); i++) {
			if (m_pop->getOffspring()[i].fitness() == 0)//连续存活n代
				count++;
		}
		if (count > m_pop->size())
			return true;
		else
			return false;
	}

	//std::vector<size_t> SPMOEA9::predictObjSpacePotential(Problem *pro) {//输出的是决策子空间的位置
	//	//在目标空间上表现为：与已收敛子空间处于互不支配地位的子空间且里面个体很少，还有局部PF的边界点周围
	//	//在搜索空间上表现为：目标空间上述点在搜索空间对应点的邻域范围内
	//	//std::vector<std::vector<size_t>> potential_spaces;
	//	std::vector<size_t> indi_index;//第一排个体索引
	//	for (size_t i = 0; i < m_pop->size(); i++) {
	//		if (m_pop->at(i).fitness() == 0)
	//			indi_index.emplace_back(i);
	//	}
	//	size_t num_obj = CAST_CONOP(pro)->numberObjectives();
	//	std::vector<size_t> obj_point_index;//第一排个体的子目标边界值索引
	//	for (size_t i = 0; i < num_obj; i++) {
	//		size_t temp = indi_index[0];
	//		for (size_t j = 0; j < indi_index.size(); j++) {
	//			if (m_pop->at(temp).objective()[i] > m_pop->at(indi_index[j]).objective()[i]) {
	//				temp = indi_index[j];
	//			}
	//		}
	//		if (obj_point_index.empty() || std::find(obj_point_index.begin(), obj_point_index.end(), temp) == obj_point_index.end())
	//			obj_point_index.emplace_back(temp);
	//	}
	//	//再加入第一排拥挤距离的索引,加入个体在单个目标上的最小距离大于平均值对应的个体，整体最小距离大于平均距离的个体
	//	std::map<size_t, std::vector<Real>> crowd_dist;//前排个体对应的拥挤距离
	//	for (size_t i = 0; i < indi_index.size(); ++i) {
	//		std::pair<size_t, std::vector<Real>> indi_dist;
	//		std::vector<Real> m_dist;
	//		for (size_t j = 0; j < num_obj + 1; ++j) {
	//			indi_dist.first = indi_index[i];
	//			Real dist = (Real)std::numeric_limits<int>::max();
	//			for (size_t k = 0; k < indi_index.size(); ++k) {
	//				if (j != num_obj) {
	//					if (k != i) {
	//						Real temp = fabs(m_pop->at(indi_index[k]).objective()[j] - m_pop->at(indi_index[i]).objective()[j]);
	//						if (temp < dist)
	//							dist = temp;
	//					}
	//				}
	//				else {
	//					if (k != i) {
	//						Real temp = euclideanDistance(m_pop->at(indi_index[k]).objective().begin(), m_pop->at(indi_index[k]).objective().end(), m_pop->at(indi_index[i]).objective().begin());
	//						if (temp < dist)
	//							dist = temp;
	//					}
	//				}
	//			}
	//			m_dist.emplace_back(dist);
	//		}
	//		indi_dist.second = m_dist;
	//		crowd_dist.insert(indi_dist);
	//	}
	//	//计算每个距离指标上的平均值
	//	std::vector<Real> average_dist;
	//	Real total_dist;
	//	for (size_t i = 0; i < num_obj + 1; ++i) {
	//		total_dist = 0;
	//		for (const auto& ind : crowd_dist) {
	//			total_dist += ind.second[i];
	//		}
	//		average_dist.emplace_back(total_dist / indi_index.size());
	//	}
	//	//得到潜力个体的索引
	//	for (size_t i = 0; i < num_obj + 1; ++i) {
	//		for (const auto& ind : crowd_dist) {
	//			if (ind.second[i] > average_dist[i]) {
	//				if (std::find(obj_point_index.begin(), obj_point_index.end(), ind.first) == obj_point_index.end())
	//					obj_point_index.emplace_back(ind.first);
	//			}
	//		}
	//	}
	//	//根据个体索引，找到搜索空间对应子空间，并为子空间设置潜力值，潜力值大小为下一次在该子空间产生子代的个数
	//	std::vector<size_t> predict_index;
	//	Real ratio = 0.9*(0.1 * m_pop->size() + 0.9 * indi_index.size() - 1) / (m_pop->size() - 1);//开发的比例
	//	for (size_t i = 0; i < obj_point_index.size(); i++) {
	//		auto position = m_pop->at(obj_point_index[i]).variable().vect();
	//		auto index = m_mo_hlc->subspaceTree().getRegionIdx(position);
	//		//size_t temp_potential = m_mo_hlc->getSubspacePotential(index);
	//		//m_mo_hlc->setSubspacePotential(index, temp_potential + m_pop->size() * ratio/obj_point_index.size());
	//		if (predict_index.empty() || std::find(predict_index.begin(), predict_index.end(), index) == predict_index.end())
	//			predict_index.emplace_back(index);
	//	}
	//	//updateVarSpacePotential();
	//	return predict_index;
	//}

	std::pair<std::vector<size_t>,std::vector<std::vector<size_t>>> SPMOEA9::predictObjSpacePotential(Problem *pro) {//输出的是决策子空间的位置
		//在目标空间上表现为：与已收敛子空间处于互不支配地位的子空间且里面个体很少，还有局部PF的边界点周围
		//在搜索空间上表现为：目标空间上述点在搜索空间对应点的邻域范围内
		std::vector<std::vector<size_t>> potential_ind_idx;//第一行为所有目标的极值，剩余的为各子目标下的稀疏点对
		std::vector<size_t> first_indi_index;//第一排个体索引
		for (size_t i = 0; i < m_pop->size(); i++) {
			if (m_pop->at(i).fitness() == 0)
				first_indi_index.emplace_back(i);
		}
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		std::vector<size_t> obj_point_index;
		//先加入前排个体中子目标最差的边界值索引
		for (size_t i = 0; i < num_obj; i++) {
			std::vector<size_t> boundary_idx;
			size_t temp = first_indi_index[0];
			for (size_t j = 0; j < first_indi_index.size(); j++) {
				if (m_pop->at(temp).objective()[i] < m_pop->at(first_indi_index[j]).objective()[i]) {
					temp = first_indi_index[j];
				}
			}
			boundary_idx.push_back(temp);

			if (obj_point_index.empty() || std::find(obj_point_index.begin(), obj_point_index.end(), temp) == obj_point_index.end()) {
				obj_point_index.emplace_back(temp);
			}
			potential_ind_idx.push_back(boundary_idx);
		}
		
		//再加入第一排稀疏个体的索引，计算距离前需要对目标空间进行归一化
		// 在每个子目标上检索距离最远的个体
		std::vector<size_t> sparse_idx;
		std::vector<std::vector<Real>> temp_first_obj;
		for (size_t i = 0; i < first_indi_index.size() ; ++i) {
			auto temp_obj = m_pop->at(first_indi_index[i]).objective();
			if (m_normalize) {
				for (size_t j = 0; j < temp_obj.size(); ++j) {
					temp_obj[j] = (temp_obj[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
				}
			}
			temp_first_obj.push_back(temp_obj);
		}
		if (first_indi_index.size() > 1) {
			for (size_t i = 0; i < num_obj; ++i) {
				//所有解按照某一维目标排序
				//先构造元组
				std::vector<size_t> sparse_pair_idx;//存放每个目标上的稀疏点对
				std::vector<std::tuple<Real, size_t>> temp_dim_obj;
				for (size_t j = 0; j < temp_first_obj.size(); ++j) {
					temp_dim_obj.emplace_back(std::make_tuple<>(temp_first_obj[j][i], first_indi_index[j]));
				}
				//按目标值递增排序
				for (size_t j = 0; j < temp_dim_obj.size(); ++j) {
					for (size_t k = j + 1; k < temp_dim_obj.size(); ++k) {
						if (std::get<0>(temp_dim_obj[j]) > std::get<0>(temp_dim_obj[k])) {
							auto temp = temp_dim_obj[j];
							temp_dim_obj[j] = temp_dim_obj[k];
							temp_dim_obj[k] = temp;
						}
					}
				}
				//求最大间隔
				Real max_dist = 0;
				size_t idx = 0;
				for (size_t j = 1; j < temp_dim_obj.size(); ++j) {
					auto temp = std::get<0>(temp_dim_obj[j]) - std::get<0>(temp_dim_obj[j-1]);
					if (temp > max_dist) {
						max_dist = temp;
						idx = j;
					}
				}
				sparse_pair_idx.push_back(std::get<1>(temp_dim_obj[idx-1]));
				sparse_pair_idx.push_back(std::get<1>(temp_dim_obj[idx]));
				//判断是否已有相同的点对加入
				if (potential_ind_idx.size() == num_obj) {
					potential_ind_idx.push_back(sparse_pair_idx);
				}
				else {
					bool flag = false;
					for (size_t j = num_obj; j < potential_ind_idx.size(); ++j) {
						auto obj = potential_ind_idx[j];
						if (obj[0] == sparse_pair_idx[0] && obj[1] == sparse_pair_idx[1] || \
							obj[0] == sparse_pair_idx[1] && obj[1] == sparse_pair_idx[0]) {
							flag = true;
							break;
						}
					}
					if (!flag) {
						potential_ind_idx.push_back(sparse_pair_idx);
					}
				}
				
				if (sparse_idx.empty()||std::find(sparse_idx.begin(), sparse_idx.end(),std::get<1>(temp_dim_obj[idx])) == sparse_idx.end()) {
					sparse_idx.push_back(std::get<1>(temp_dim_obj[idx]));
				}
				if (sparse_idx.empty() || std::find(sparse_idx.begin(), sparse_idx.end(), std::get<1>(temp_dim_obj[idx-1])) == sparse_idx.end()) {
					sparse_idx.push_back(std::get<1>(temp_dim_obj[idx-1]));
				}
			}
		}
		for (size_t i = 0; i < sparse_idx.size(); ++i) {
			if (std::find(obj_point_index.begin(), obj_point_index.end(), sparse_idx[i]) == obj_point_index.end()) {
				obj_point_index.push_back(sparse_idx[i]);
			}
		}

		//最后根据个体索引得到搜索子空间索引
		std::vector<size_t> var_space_idx;
		for (size_t i = 0; i < obj_point_index.size(); ++i) {
			auto position = m_pop->at(obj_point_index[i]).variable().vect();
			auto index = m_mo_hlc->subspaceTree().getRegionIdx(position);
			if (var_space_idx.empty() || std::find(var_space_idx.begin(), var_space_idx.end(), index) == var_space_idx.end()) {
				var_space_idx.push_back(index);
			}
		}
		std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> output_pair;
		output_pair.first = var_space_idx;
		output_pair.second = potential_ind_idx;
		return output_pair;
	}

	//void SPMOEA9::updateVarSpacePotential() {
	//	std::vector<size_t> sub_potential;
	//	size_t sum = 0;
	//	//for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
	//	//	sub_potential.push_back(m_mo_hlc->getSubspaceInfo(i).m_sub_potential);
	//	//	sum += sub_potential.back();
	//	//}
	//	//Real ratio = 0.9;//开发比例
	//	//for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
	//	//	m_mo_hlc->setSubspacePotential(i,std::floor(ratio*m_pop->size()*sub_potential[i]/sum));
	//	//}
	//}

	/*std::vector<size_t> SPMOEA9::predictVarSpacePotential(Problem *pro) {
		std::vector<size_t> space_potential;

		return space_potential;
	}*/

	void SPMOEA9::envirSelection(Problem *pro, Random *rnd) {//是否趋近收敛有不同的环境选择策略
		//m_pop->envirSelection(pro, rnd, b);
		//基于子空间跨度的多参考选择方法
		if (m_subspace_select) {
			auto select_pop = selectIndi(m_pop->getOffspring(), m_pop->size(), pro, rnd);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_pop->at(i) = m_pop->getOffspring()[select_pop[i]];
				m_pop->at(i).setCounter(m_pop->at(i).surviveAge() + 1);
			}
		}
		else {
			m_pop->envirSelection(pro, rnd, false);
		}
	}

	void SPMOEA9_pop::envirSelection(Problem *pro, Random *rnd, bool b) {
		if (b) {//趋近收敛时的环境选择基于子空间
			//基于目标子空间进行选择，选择依据为：前沿子空间的均匀性和支配空间的收敛性
			//先对父代和子代进行非支配排序
			/*for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				m_offspring[i + this->m_individuals.size()] = *this->m_individuals[i];
			}*/
			nondominatedSorting(m_offspring);
			/*auto select_pop = selectIndi(m_offspring, this->m_individuals.size(),pro);
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				*this->m_individuals[i] = m_offspring[select_pop[i]];
				this->m_individuals[i]->setCounter(this->m_individuals[i]->surviveAge() + 1);
			}*/
		}
		else {
			survivorSelection(*this, m_offspring);
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				this->m_individuals[i]->setCounter(this->m_individuals[i]->surviveAge() + 1);
			}
		}
	}

	//std::vector<Real> SPMOEA9_pop::optimaPredict(const std::vector<std::vector<Real>>& points) {//预测得到的是坐标,points为一个点序列
	//	std::vector<Real> m_coordinate;
	//	std::vector<std::vector<Real>> m_2dim_point;//每一个二维的预测点
	//	size_t num = points[0].size() - 1;//num个二维预测
	//	for (size_t i = 0; i < num; ++i) {
	//		std::vector<Real> m_length(points.size() - 1);//一维预测数据
	//		std::vector<Real> m_angle(points.size() - 1);
	//		std::vector<Real> m_begin = { points[0][0],points[0][i + 1] };
	//		for (size_t j = 0; j < points.size() - 1; ++j) {
	//			m_length[j] = sqrt(pow(points[j + 1][0] - points[j][0], 2) + pow(points[j + 1][i + 1] - points[j][i + 1], 2));
	//			m_angle[j] = asin((points[j + 1][i + 1] - points[j][i + 1]) / m_length[j]);
	//		}
	//		//采用最小二乘拟合得到预测点位置
	//		auto temp1 = least_square_estimation(m_length);
	//		auto temp2 = least_square_estimation(m_angle);
	//		std::vector<Real> predict_point;
	//		predict_point.emplace_back(points.back()[0] + temp1 * cos(temp2));
	//		predict_point.emplace_back(points.back()[i + 1] + temp1 * sin(temp2));
	//		m_2dim_point.emplace_back(predict_point);
	//	}
	//	for (size_t i = 0; i < m_2dim_point.size() - 1; ++i) {
	//		m_2dim_point[i + 1][1] = m_2dim_point[i + 1][1] * m_2dim_point[0][0] / m_2dim_point[i + 1][0];
	//	}
	//	for (size_t i = 0; i < m_2dim_point.size() + 1; ++i) {
	//		if (i == 0)
	//			m_coordinate.emplace_back(m_2dim_point[0][0]);
	//		else if (i == 1)
	//			m_coordinate.emplace_back(m_2dim_point[0][0]);
	//		else
	//			m_coordinate.emplace_back(m_2dim_point[i - 1][1]);
	//	}
	//	return m_coordinate;
	//}

	//std::vector<std::vector<std::vector<Real>>> SPMOEA9_pop::indi_time_series(const std::vector<std::vector<Solution<>>>& p) {//对memeory中个体直接找距离最小值形成序列
	//	std::vector<std::vector<std::vector<Real>>> all_series;
	//	for (size_t i = 0; i < p[0].size(); ++i) {
	//		std::vector<std::vector<Real>> one_series;
	//		for (size_t j = 0; j < p.size(); ++j) {
	//			if (j == 0) {
	//				one_series.emplace_back(p[0][i].variable().vect());
	//			}
	//			else {
	//				Real dist = euclideanDistance(one_series.back().begin(), one_series.back().end(), p[j][0].variable().vect().begin());
	//				size_t count = 0;
	//				for (size_t k = 1; k < p[j].size(); ++k) {
	//					Real temp = euclideanDistance(one_series.back().begin(), one_series.back().end(), p[j][k].variable().vect().begin());
	//					if (dist > temp)
	//						count = k;
	//				}
	//				one_series.emplace_back(p[j][count].variable().vect());
	//			}
	//		}
	//		all_series.emplace_back(one_series);
	//	}
	//	return all_series;
	//}

	std::vector<size_t> SPMOEA9::selectIndi(const Population<Solution<>>& pop, size_t select_num,Problem *pro, Random *rnd) {
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		auto pop_att = popAttach(pop);
		std::vector<size_t> select_index;//已经选择个体的索引
		std::vector<size_t> behind_space_index;//后排子空间索引
		std::vector<size_t> front_space_index;//前排子空间索引
		size_t front_space_num = 0;//前沿子空间内的个体总数
		std::vector<size_t> first_ind_index;//先得到第一排个体的索引
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		for (const auto& ind : pop_att[0][0]) {
			front_space_index.emplace_back(ind.first);
			front_space_num += ind.second.size();
		}
		for (const auto& ind : pop_att[0][1]) {
			behind_space_index.emplace_back(ind.first);
		}
		std::map<size_t, std::vector<size_t>> map_first_indi;//得到前排个体分布在哪些子空间
		for (const auto& sp : pop_att[0][0]) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = sp.first;
			for (const auto& in : sp.second) {
				if (pop[in].fitness() == 0) {
					temp.second.push_back(in);
				}
			}
			map_first_indi.insert(temp);
		}
		std::map<size_t, std::vector<size_t>> map_selected_indi;//前沿子空间已经选的解，顺序索引
		std::map<size_t, std::vector<size_t>> map_no_selected_indi;//前沿子空间已经选的解，顺序索引
		//根据前排子空间内部的个体总数进行个体选择，选择时考虑个体年龄
		if (first_ind_index.size() > m_pop->size()) {//前沿个体充足时，进行子空间多参考点选择
			//先得到当前子目标的极值解的索引
			std::vector<size_t> sub_opt_idx;
			for (size_t j = 0; j < num_obj; ++j) {
				std::vector<Real> obj_value;
				for (size_t k = 0; k < first_ind_index.size(); ++k) {
					obj_value.push_back(pop[first_ind_index[k]].objective()[j]);
				}
				size_t min_idx = min_element(obj_value.begin(), obj_value.end()) - obj_value.begin();
				sub_opt_idx.push_back(min_idx);
			}
		    //再选择子空间内某个子目标的极值点
			for (const auto& i : map_first_indi) {
				std::pair<size_t, std::vector<size_t>> temp1;
				std::pair<size_t, std::vector<size_t>> temp2;
				temp1.first = i.first;
				temp2.first = i.first;
				auto ind_idx = i.second;
				if (ind_idx.size() == 1) {
					select_index.push_back(ind_idx[0]);
					temp1.second.push_back(ind_idx[0]);
				}
				else {//比较子空间内点子目标值
					std::vector<Real> obj_value;
					for (size_t k = 0; k < ind_idx.size(); ++k) {
						obj_value.push_back(pop[ind_idx[k]].objective()[0]);
					}
					size_t min_idx = min_element(obj_value.begin(), obj_value.end()) - obj_value.begin();
					select_index.push_back(ind_idx[min_idx]);
					temp1.second.push_back(min_idx);
				}
				std::vector<size_t> no_select_idx;
				for (size_t j = 0; j < map_first_indi[i.first].size(); ++j) {
					if (std::find(temp1.second.begin(), temp1.second.end(),j) == temp1.second.end()) {
						no_select_idx.push_back(j);
					}
				}
				temp2.second = no_select_idx;
				map_selected_indi.insert(temp1);
				map_selected_indi.insert(temp2);
			}
			std::map<size_t, std::vector<std::vector<Real>>> map_subspace_dist_matrix;//子空间未选点与其邻域已选点的距离矩阵
			std::map<size_t, std::vector<size_t>> map_neigh_select_ind;//邻域已选个体索引,真实索引
			std::map<size_t, std::vector<size_t>> map_neigh_inx;//邻域索引
			for (const auto& i : map_first_indi) {
				std::pair<size_t, std::vector<std::vector<Real>>> temp;
				std::pair<size_t, std::vector<size_t>> neigh_select_ind;
				std::pair<size_t, std::vector<size_t>> neigh_inx;
				neigh_select_ind.first = i.first;
				temp.first = i.first;
				neigh_inx.first = i.first;
				auto idx = i.second;
				//先找出邻域子空间索引，并得到所有邻域已选个体索引
				std::list<size_t> neighbors;
				std::vector<size_t> front_nei_select_idx;//邻域已选个体，实际索引
				std::vector<size_t> front_nei_idx;
				m_mo_hlc->getObjspaceTree().findNeighbor(i.first, neighbors);
				for (const auto& ss : map_first_indi) {
					if (find(neighbors.begin(), neighbors.end(), ss.first) != neighbors.end()) {
						front_nei_idx.push_back(ss.first);
						//front_nei_select_idx.push_back(ss.first);
						auto inx = ss.second;
						auto ind_selected = map_selected_indi[ss.first];
						for (size_t j = 0; j < ind_selected.size(); ++j) {
							front_nei_select_idx.push_back(ss.second[ind_selected[j]]);
						}
					}
				}
				neigh_inx.second = front_nei_idx;
				map_neigh_inx.insert(neigh_inx);
				//子空间未选个体索引
				std::vector<size_t> no_select_inx;
				auto no_se = map_no_selected_indi[i.first];
				for (size_t j = 0; j < no_se.size(); ++j) {
					no_select_inx.push_back(i.second[no_se[j]]);
				}
				//未选点与已选点的距离矩阵，如果距离矩阵为空怎么办？
				std::vector<std::vector<Real>> dist_matrix;
				for (size_t j = 0; j < no_select_inx.size(); ++j) {
					std::vector<Real> ind_dist;
					for (size_t k = 0; k < front_nei_select_idx.size(); ++k) {
						ind_dist.push_back(euclideanDistance(pop[no_select_inx[j]].objective().begin(), pop[no_select_inx[j]].objective().end(), pop[front_nei_select_idx[k]].objective().begin()));
					}
					dist_matrix.push_back(ind_dist);
				}
				temp.second = dist_matrix;
				neigh_select_ind.second = front_nei_select_idx;
				map_subspace_dist_matrix.insert(temp);
				map_neigh_select_ind.insert(neigh_select_ind);
			}
			//然后根据每个子空间内个体最小距离的统计指标选择要选的子空间
			size_t select_subspace_idx=0;
			while (select_index.size()<select_num) {
				//选择子空间内未选个体最小距离最大的子空间，从里面采用最小距离最大规则选出一个解
				std::vector<std::pair<size_t,Real>> sub_min_max;
				std::map<size_t,size_t> sub_min_max_inx;//每个子空间最大值的索引
				for (const auto& i : map_subspace_dist_matrix) {
					auto matri = i.second;
					Real max_dist = 0.;
					size_t inx = 0;
					if (matri.size() != 0) {
						std::vector<size_t> temp_min;
						for (size_t j = 0; j < i.second.size(); ++j) {
							auto min_temp_dist = *min_element(i.second[j].begin(), i.second[j].end());
							temp_min.push_back(min_temp_dist);
						}
						max_dist = *max_element(temp_min.begin(), temp_min.end());
						inx = max_element(temp_min.begin(), temp_min.end()) - temp_min.begin();
					}
					else {
						size_t a = 1;
					}
					std::pair<size_t, Real> temp1;
					std::pair<size_t, size_t> temp2;
					temp1.first = i.first;
					temp1.second = max_dist;
					temp2.first = i.first;
					temp2.second = inx;
					sub_min_max.push_back(temp1);
					sub_min_max_inx.insert(temp2);
				}
				//选出子空间
				size_t subspace_idx;
				std::vector<std::pair<size_t, Real>> subspace_dist= sub_min_max;
				for (size_t j = 0; j < sub_min_max.size()-1; ++j) {
					for (size_t k = j+1; k < sub_min_max.size(); ++k) {
						if (subspace_dist[j].second < subspace_dist[k].second) {
							std::pair<size_t, Real> temp_pair;
							temp_pair.first = subspace_dist[k].first;
							temp_pair.second = subspace_dist[k].second;
							subspace_dist[k].first= subspace_dist[j].first;
							subspace_dist[k].second = subspace_dist[j].second;
							subspace_dist[j].first = temp_pair.first;
							subspace_dist[j].second = temp_pair.second;
						}
					}
				}
				for (size_t i = 0; i < subspace_dist.size(); ++i) {
					if (map_selected_indi[subspace_dist[i].first].size() < map_first_indi[subspace_dist[i].first].size()) {
						subspace_idx = subspace_dist[i].first;
						break;
					}
				}
				size_t max_idx = sub_min_max_inx[subspace_idx];//指示未选点的索引
				auto temp_idx = map_no_selected_indi[subspace_idx][max_idx];
				//所选点的索引
				size_t sele_idx = map_first_indi[subspace_idx][temp_idx];//未选点的实际索引
				//更新邻域子空间所选点和距离矩阵
				auto nei_idx = map_neigh_inx[subspace_idx];
				for (const auto& i : map_first_indi) {
					if (std::find(nei_idx.begin(), nei_idx.end(), i.first) != nei_idx.end()) {
						if (i.first == subspace_idx) {
							//去除选中点所在行
							map_subspace_dist_matrix[i.first].erase(map_subspace_dist_matrix[i.first].begin() + max_idx);
							map_no_selected_indi[i.first].erase(map_no_selected_indi[i.first].begin() + max_idx);
							map_selected_indi[i.first].push_back(temp_idx);
						}
						//添加已选点的列
						std::vector<size_t> no_sele_ind;
						for (size_t j = 0; j < map_no_selected_indi[i.first].size(); ++j) {
							no_sele_ind.push_back(i.second[map_no_selected_indi[i.first][j]]);
						}
						for (size_t j = 0; j < map_subspace_dist_matrix[i.first].size(); ++j) {
							map_subspace_dist_matrix[i.first][j].push_back(euclideanDistance(pop[no_sele_ind[j]].objective().begin(), pop[no_sele_ind[j]].objective().end(), pop[sele_idx].objective().begin()));
						}
						map_neigh_select_ind[i.first].push_back(sele_idx);
					}
				}
				select_index.push_back(sele_idx);
			}
		}
		else if (front_space_num > m_pop->size()) {//前沿个体不足时，前沿子空间个体充足时
			bool flag = 1;
			if (flag) {//目标空间选择
				//先选前沿个体
				select_index = first_ind_index;
				//再依次选择前排子空间中的其他个体
				auto front_subspace_ind = pop_att[0][0];
				for (const auto& i : front_subspace_ind) {
					std::pair<size_t, std::vector<size_t>> temp;
					temp.first = i.first;
					std::vector<size_t> f_ind_idx;
					for (size_t j = 0; j < i.second.size(); ++j) {
						if (find(first_ind_index.begin(), first_ind_index.end(),i.second[j]) != first_ind_index.end()) {
							f_ind_idx.push_back(i.second[j]);
						}
					}
					temp.second = f_ind_idx;
					map_selected_indi.insert(temp);//前沿子空间已经选的解
				}
				while (select_index.size()<select_num) {
					for (const auto& i : front_subspace_ind) {
						if (map_selected_indi[i.first].size() < i.second.size()) {
							//从剩余个体中选择一个
							for (size_t j = 0; j < i.second.size(); ++j) {
								if (find(map_selected_indi[i.first].begin(), map_selected_indi[i.first].end(), i.second[j]) == map_selected_indi[i.first].end()) {
									select_index.push_back(i.second[j]);
									map_selected_indi[i.first].push_back(i.second[j]);
									break;
								}
							}
						}
						if (select_index.size() >= select_num) {
							break;
						}
					}
				}
			}
			else {//解空间选择，根据解空间分布选择，或者根据子空间的潜力值选择
				//将前沿子空间的点分成两部分，前沿解和非前沿解
				
				//将所有解映射到搜索空间，找出离前沿解最远的个体
			}
		}
		else {//前排个体数不足时，前排全选，后面的轮流选
			//先加入前排子空间个体的索引
			for (const auto& ind : pop_att[0][0]) {
				for (const auto& p : ind.second) {
					select_index.push_back(p);
				}
			}
			//再在后排子空间轮流选择
			std::vector<size_t> space_select_num(pop_att[0][1].size(), 0);//每个后排子空间已选择的个数
			std::vector<std::vector<size_t>> behind_spaces;
			for (const auto& sp : pop_att[0][1]) {
				behind_spaces.push_back(sp.second);
			}
			while (select_index.size() < select_num) {
				for (size_t i = 0; i < behind_spaces.size(); ++i) {
					if (space_select_num[i] < behind_spaces.size()) {
						size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
						if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
							select_index.push_back(behind_spaces[i][select_idx]);
							space_select_num[i]++;
						}
					}
					if (select_index.size() >= select_num) {
						break;
					}
				}
			}
		}
		//if (front_space_num > m_pop->size()) {//前排个体充足时，进行子空间多参考点选择
		//	//计算前排子空间跨度
		//	std::map<size_t, size_t> space_select_num;
		//	std::vector<Real> subspace_span;
		//	for (const auto& subspace : pop_att[0][0]) {
		//		std::vector<Real> ideal_point;
		//		std::vector<Real> nadir_point;
		//		auto sub_boundary = m_obj_tree->getBox(subspace.first);
		//		for (size_t j = 0; j < sub_boundary.size(); ++j) {
		//			ideal_point.push_back(sub_boundary[j].first);
		//			nadir_point.push_back(sub_boundary[j].second);
		//		}
		//		std::vector<size_t> ind_index = subspace.second;
		//		std::vector<std::vector<Real>> sub_points;
		//		for (size_t q = 0; q < ind_index.size(); ++q) {
		//			auto temp = pop[ind_index[q]].objective();
		//			if (m_normalize) {
		//				for (size_t j = 0; j < temp.size(); ++j) {
		//					temp[j] = (temp[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
		//				}
		//			}
		//			sub_points.push_back(temp);
		//		}
		//		subspace_span.push_back(calSubspaceSpan(ideal_point, nadir_point, sub_points));
		//	}
		//	//为前排子空间分配个体数
		//	std::vector<size_t> subspace_ind_num;
		//	Real sum_span = 0.;
		//	for (size_t p = 0; p < subspace_span.size(); ++p) {
		//		sum_span += subspace_span[p];
		//	}
		//	for (size_t p = 0; p < subspace_span.size(); ++p) {
		//		subspace_ind_num.push_back(1 + std::ceil(m_pop->size() * subspace_span[p] / sum_span));
		//	}
		//	//从前排子空间选择个体，满足优先级，即先子目标极值，再前排个体，最后其他个体
		//	size_t index = 0;
		//	m_over_space_index.clear();
		//	for (const auto& subspace : map_first_indi) {
		//		std::vector<Real> anchor_point;
		//		auto sub_boundary = m_obj_tree->getBox(subspace.first);
		//		for (size_t p = 0; p < sub_boundary.size(); ++p) {
		//			anchor_point.push_back(sub_boundary[p].first);
		//		}
		//		std::vector<size_t> ind_index = subspace.second;
		//		/*std::vector<std::vector<Real>> sub_points;
		//		for (size_t p = 0; p < ind_index.size(); ++p) {
		//			sub_points.push_back(pop[ind_index[p]].objective());
		//		}*/
		//		auto select_ind=selectIndSubspace(anchor_point, ind_index, subspace_ind_num[index]);
		//		std::vector<size_t> select_ind_index = std::get<0>(select_ind);
		//		bool flag= std::get<1>(select_ind);
		//		if (flag) {
		//			m_over_space_index.push_back(subspace.first);
		//		}
		//		for (size_t p = 0; p < select_ind_index.size(); ++p) {
		//			select_index.push_back(select_ind_index[p]);
		//		}
		//		index++;
		//		space_select_num.insert(std::make_pair<>(subspace.first,select_ind_index.size()));
		//	}
		//	//当选择个体超出时，进行拥挤淘汰，1216行越界
		//	while (select_index.size() > m_pop->size())
		//	{
		//		//先找出子目标极值个体
		//		std::vector<size_t> sub_optima_inx;
		//		for (size_t p = 0; p < num_obj; ++p) {
		//			size_t inx = 0;
		//			Real min_v = 1. * 10e14;
		//			for (size_t q = 0; q < select_index.size(); ++q) {
		//				if (pop[select_index[q]].objective()[p] < min_v) {
		//					min_v = pop[select_index[q]].objective()[p];
		//					inx = q;
		//				}
		//			}
		//			if (sub_optima_inx.empty() || std::find(select_index.begin(), select_index.end(), inx) == select_index.end()) {
		//				sub_optima_inx.push_back(select_index[inx]);
		//			}
		//		}
		//		//再计算个体的拥挤度，只计算与邻域个体的距离，有问题
		//		std::vector<Real> density;
		//		for (size_t p = 0; p < select_index.size(); ++p) {
		//			Real min_dist = 1.*10e14;
		//			if (std::find(sub_optima_inx.begin(), sub_optima_inx.end(), select_index[p]) == sub_optima_inx.end()) {
		//				//索引p所在子空间的的邻域中选中的个体之间的距离
		//				std::vector<size_t> compare_idx;
		//				auto temp_obj = pop[select_index[p]].objective();
		//				if (m_normalize) {
		//					for (size_t m = 0; m < temp_obj.size(); ++m) {
		//						temp_obj[m] = (temp_obj[m] - m_front_pop_range[m].first) / (m_front_pop_range[m].second - m_front_pop_range[m].first);
		//					}
		//				}
		//				size_t sp_idx = m_obj_tree->getRegionIdx(temp_obj);
		//				std::list<size_t> neighbors;
		//				m_obj_tree->findNeighbor(sp_idx,neighbors);
		//				neighbors.push_back(sp_idx);
		//				for (const auto ne : neighbors) {
		//					if (std::find(front_space_index.begin(), front_space_index.end(), ne) != front_space_index.end()) {
		//						for (const auto& ind : pop_att[0][0][ne]) {
		//							if (std::find(select_index.begin(), select_index.end(), ind) != select_index.end()) {
		//								compare_idx.push_back(ind);
		//							}
		//						}
		//					}
		//				}
		//				for (size_t k = 0; k < compare_idx.size(); ++k) {
		//					if (p != k) {
		//						auto p1 = pop[select_index[p]].objective();
		//						auto p2 = pop[compare_idx[k]].objective();
		//						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
		//						if (dist < min_dist) {
		//							min_dist = dist;
		//						}
		//					}
		//				}
		//			}
		//			density.push_back(min_dist);
		//		}
		//		//再去掉最拥挤的
		//		size_t inx = 0;
		//		Real dist = 1.*10e14;
		//		for (size_t p = 0; p < density.size(); ++p) {
		//			if (density[p] < dist) {
		//				dist = density[p];
		//				inx = p;
		//			}
		//		}
		//		select_index.erase(select_index.begin() + inx);
		//	}

		//	//当选择总数不足时,依次从前排子空间中轮流选择剩余的个体，标记选完的子空间
		//	std::map<size_t, std::vector<size_t>> all_ind_map= pop_att[0][0];
		//	
		//	while (select_index.size() < m_pop->size()) {
		//		for (const auto& sp : all_ind_map) {
		//			if (space_select_num[sp.first] < sp.second.size()) {//未选完的子空间
		//				while (select_index.size() < m_pop->size()) {
		//					for (const auto& ind : sp.second) {
		//						if (std::find(select_index.begin(), select_index.end(), ind) == select_index.end()) {
		//							select_index.push_back(ind);
		//							space_select_num[sp.first]++;
		//							break;
		//						}
		//					}
		//					break;
		//				}
		//			}
		//		}
		//	}
		//}
		//else {//前排个体数不足时，前排全选，后面的轮流选
		//	//先加入前排子空间个体的索引
		//	for (const auto &ind :pop_att[0][0]) {
		//		for (const auto& p : ind.second) {
		//			select_index.push_back(p);
		//		}
		//	}
		//	//再在后排子空间轮流选择
		//	std::vector<size_t> space_select_num(pop_att[0][1].size(),0);//每个后排子空间已选择的个数
		//	std::vector<std::vector<size_t>> behind_spaces;
		//	for (const auto& sp : pop_att[0][1]){
		//		behind_spaces.push_back(sp.second);
		//	}
		//	while (select_index.size() < m_pop->size()) {
		//		for (size_t i = 0; i < behind_spaces.size();++i) {
		//			if (space_select_num[i] < behind_spaces.size()) {
		//				size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
		//				if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
		//					select_index.push_back(behind_spaces[i][select_idx]);
		//					space_select_num[i]++;
		//				}
		//			}
		//			if (select_index.size() >= m_pop->size()) {
		//				break;
		//			}
		//		}
		//	}
		//}
		return select_index;
	}

	//Real SPMOEA9::calSubspaceSpan(const std::vector<Real>& ideal_point, const std::vector<Real>& nadir_point, const std::vector<std::vector<Real>>& points) {
	//	//先判断子空间流形的凹凸，决定选择哪一个为参考点
	//	//先找边界点
	//	if (points.size() == 1) {
	//		return 0.00001;
	//	}
	//	else {
	//		std::vector<std::vector<Real>> boundary_points;
	//		size_t num_obj = ideal_point.size();
	//		for (size_t i = 0; i < num_obj; ++i) {
	//			if (i + 1 <= points.size()) {
	//				size_t inx = 0;
	//				Real min_v = 1. * 10e14;
	//				for (size_t j = 0; j < points.size(); ++j) {
	//					if (points[j][i] < min_v) {
	//						min_v = points[j][i];
	//						inx = j;
	//					}
	//				}
	//				boundary_points.push_back(points[inx]);
	//			}
	//		}

	//		bool if_convex = false;
	//		//记录边界点到理想点的最小值
	//		Real min_v = 1. * 10e14;
	//		for (const auto& p : boundary_points) {
	//			auto dist = euclideanDistance(ideal_point.begin(), ideal_point.end(), p.begin());
	//			if (dist < min_v) {
	//				min_v = dist;
	//			}
	//		}
	//		//再看其他点的距离有没有小于最小值的
	//		size_t count = 0;//记数满足条件数
	//		for (const auto& p : points) {
	//			auto dist = euclideanDistance(ideal_point.begin(), ideal_point.end(), p.begin());
	//			if (dist < min_v) {
	//				count++;
	//			}
	//		}
	//		if (count > 1) {
	//			if_convex = true;
	//		}
	//		std::vector<Real> ref_point;
	//		if (if_convex) {
	//			ref_point = nadir_point;
	//		}
	//		else {
	//			ref_point = ideal_point;
	//		}
	//		auto span = manifoldSpan(ref_point, points);
	//		if (span.second < 0) {
	//			size_t a = 1;
	//		}
	//		return span.second;
	//	}
	//}

	//std::tuple<std::vector<size_t>,bool> SPMOEA9::selectIndSubspace(const std::vector<Real>& ref_point, const std::vector<size_t>& ind_index, int num) {
	//	//满足优先级，即先子目标极值，再前排个体，最后其他个体
	//    //对没选够的子空间记录
	//	std::vector<size_t> select_ind;
	//	bool complete_flag = false;
	//	size_t count = 0;
	//	//先得到前排个体的数量
	//	size_t front_ind_num = 0;
	//	std::vector<size_t> front_ind_idx;//前排个体的索引
	//	for (size_t i = 0; i < ind_index.size(); ++i) {
	//		if (m_pop->getOffspring()[ind_index[i]].fitness() == 0) {
	//			front_ind_num++;
	//			front_ind_idx.push_back(ind_index[i]);
	//		}
	//	}
	//	
	//	if (front_ind_num <= num) {
	//		//先加入前排个体，再加入子空间内其他个体，选择方式为：
	//		select_ind = front_ind_idx;
	//		while (ind_index.size() > select_ind.size() && select_ind.size() < num) {
	//			size_t inx = 0;
	//			Real max_ind_angle = -1. * 10e14;
	//			std::vector<Real> ind_min_angle;//选用角度还是弧长还是跨度计算
	//			for (size_t i = 0; i < ind_index.size(); ++i) {
	//				Real min_angle = 1. * 10e14;
	//				for (size_t j = 0; j < select_ind.size(); ++j) {
	//					auto p1 = m_pop->getOffspring()[ind_index[i]].objective();
	//					auto p2 = m_pop->getOffspring()[select_ind[j]].objective();
	//					if (m_normalize) {
	//						for (size_t k = 0; k < p1.size(); ++k) {
	//							p1[k] = (p1[k] - m_front_pop_range[k].first) / (m_front_pop_range[k].second - m_front_pop_range[k].first);
	//							p2[k] = (p2[k] - m_front_pop_range[k].first) / (m_front_pop_range[k].second - m_front_pop_range[k].first);
	//						}
	//					}
	//					auto angle_v = vectorAngle(ref_point, p1, p2);
	//					if (angle_v < min_angle) {
	//						min_angle = angle_v;
	//					}
	//				}
	//				ind_min_angle.push_back(min_angle);
	//			}
	//			for (size_t i = 0; i < ind_min_angle.size(); ++i) {
	//				if (ind_min_angle[i] > max_ind_angle) {
	//					max_ind_angle = ind_min_angle[i];
	//					inx = ind_index[i];
	//				}
	//			}
	//			if (select_ind.empty() || std::find(select_ind.begin(), select_ind.end(), inx) == select_ind.end()) {
	//				select_ind.push_back(inx);
	//			}
	//		}
	//		if (select_ind.size() < num) {
	//			complete_flag = true;
	//		}
	//	}
	//	else {
	//		std::vector<size_t> select_flag(front_ind_idx.size(), 0);//标记选择
	//		//先加入前沿点的边界点
	//		size_t num_obj = ref_point.size();
	//		for (size_t i = 0; i < num_obj; ++i) {
	//			if (count < num) {
	//				size_t inx = 0;
	//				Real max_v = -1. * 10e14;
	//				for (size_t j = 0; j < front_ind_idx.size(); ++j) {
	//					if (m_pop->getOffspring()[front_ind_idx[j]].objective()[i] > max_v) {
	//						max_v = m_pop->getOffspring()[front_ind_idx[j]].objective()[i];
	//						inx = j;
	//					}
	//				}
	//				if (select_ind.empty() || std::find(select_ind.begin(), select_ind.end(), front_ind_idx[inx]) == select_ind.end()) {
	//					select_ind.push_back(front_ind_idx[inx]);
	//					select_flag[inx] = 1;
	//					count++;
	//				}
	//			}
	//		}
	//		//再基于已选择点选择中间点
	//		std::vector<Real> ideal_point;
	//		std::vector<Real> nadir_point;
	//		size_t box_num=0;
	//		auto obj = m_pop->getOffspring()[select_ind[0]].objective();
	//		if (m_normalize) {
	//			for (size_t i = 0; i < obj.size(); ++i) {
	//				obj[i] = (obj[i] -m_front_pop_range[i].first) / (m_front_pop_range[i].second-m_front_pop_range[i].first);
	//			}
	//		}
	//		box_num = m_mo_hlc->subspaceTree().getRegionIdx(obj);
	//		auto box = m_mo_hlc->subspaceTree().getBox(box_num);
	//		for (size_t i = 0; i < box.size(); ++i) {
	//			ideal_point.push_back(box[i].first);
	//			nadir_point.push_back(box[i].second);
	//		}
	//		while (count < num) {
	//			size_t inx = 0;
	//			Real max_ind_span = -1. * 10e14;
	//			std::vector<Real> ind_min_span;//选用角度还是弧长还是跨度计算
	//			for (size_t i = 0; i < front_ind_idx.size(); ++i) {
	//				if (select_flag[i] == 0) {
	//					Real min_span = 1. * 10e14;
	//					for (size_t j = 0; j < select_ind.size(); ++j) {
	//						auto p1 = m_pop->getOffspring()[front_ind_idx[i]].objective();
	//						auto p2 = m_pop->getOffspring()[select_ind[j]].objective();
	//						if (m_normalize) {
	//							for (size_t k = 0; k < p1.size(); ++k) {
	//								p1[k] = (p1[k] - m_front_pop_range[k].first) / (m_front_pop_range[k].second - m_front_pop_range[k].first);
	//								p2[k] = (p2[k] - m_front_pop_range[k].first) / (m_front_pop_range[k].second - m_front_pop_range[k].first);
	//							}
	//						}
	//						//auto span_v = vectorAngle(ref_point, p1, p2);
	//						std::vector<std::vector<Real>> sub_point;
	//						sub_point.push_back(p1);
	//						sub_point.push_back(p2);
	//						auto span_v=calSubspaceSpan(ideal_point, nadir_point, sub_point);
	//						//auto span_v = euclideanDistance(p1.begin(),p1.end(),p2.begin());
	//						if (span_v < min_span) {
	//							min_span = span_v;
	//						}
	//					}
	//					ind_min_span.push_back(min_span);
	//				}
	//				else {
	//					ind_min_span.push_back(0.);
	//				}
	//			}
	//			Real sum = 0.;
	//			for (auto& an : ind_min_span) {
	//				sum += an;
	//			}
	//			if (sum == 0.) {
	//				for (size_t i = 0; i < front_ind_idx.size(); ++i) {
	//					if (select_ind.empty() || std::find(select_ind.begin(), select_ind.end(), front_ind_idx[i]) == select_ind.end()) {
	//						select_ind.push_back(front_ind_idx[i]);
	//						select_flag[i] = 1;
	//						count++;
	//						complete_flag = true;
	//						break;
	//					}
	//				}
	//			}
	//			else {
	//				size_t idx = 0;
	//				for (size_t i = 0; i < ind_min_span.size(); ++i) {
	//					if (ind_min_span[i] > max_ind_span) {
	//						max_ind_span = ind_min_span[i];
	//						inx = front_ind_idx[i];
	//						idx = i;
	//					}
	//				}
	//				if (select_ind.empty() || std::find(select_ind.begin(), select_ind.end(), inx) == select_ind.end()) {
	//					select_ind.push_back(inx);
	//					select_flag[idx] = 1;
	//					count++;
	//				}
	//			}
	//			
	//		}
	//	}
	//	return std::make_tuple<>(select_ind,complete_flag);
	//}

	space_attach SPMOEA9::spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b) {
		space_attach space_all_indi;
		//先得到第一排个体的索引
		std::vector<size_t> first_ind_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		//然后得到每个个体在哪个子空间
		std::vector<size_t> space_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			std::vector<Real> temp;
			if (b) {
				temp = pop[i].objective();
				if (m_normalize) {
					for (size_t j = 0; j < temp.size(); ++j) {
						temp[j] = (temp[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
					}
				}
			}
			else {
				temp = pop[i].variable().vect();
			}
			for (size_t j = 0; j < temp.size(); ++j) {
				if (temp[j] > 1) {
					size_t a = 1;//个体映射到子空间会越界，但是索引为最近的一个
				}
			}
			space_index.emplace_back(tree.getRegionIdx(temp));
		}
		//得到前排子空间索引
		std::vector<size_t> nondominate_space;
		for (size_t i = 0; i < first_ind_index.size(); ++i) {
			size_t temp1 = space_index[first_ind_index[i]];
			if (nondominate_space.empty())
				nondominate_space.emplace_back(temp1);
			else if (std::find(nondominate_space.begin(), nondominate_space.end(), temp1) == nondominate_space.end()) {
				nondominate_space.emplace_back(temp1);
			}
		}
		//前排子空间内的个体总数
		size_t front_indi_num = 0;
		//得到前排子空间内包含的所有个体
		std::map<size_t, std::vector<size_t>> front_indi;//key为子空间索引，value为个体索引
		std::vector<size_t> redi_space_index = space_index;
		for (size_t i = 0; i < nondominate_space.size(); ++i) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = nondominate_space[i];
			for (size_t j = 0; j < space_index.size(); ++j) {
				if (space_index[j] == temp.first) {
					temp.second.emplace_back(j);
					++front_indi_num;
				}
			}
			front_indi.insert(temp);
		}
		space_all_indi.emplace_back(front_indi);
		
		//非前排子空间索引
		std::vector<size_t> dominated_space;
		//得到非前排子空间内包含的个体
		//std::vector<std::pair<size_t, std::vector<size_t>>> behind_indi;
		std::map<size_t, std::vector<size_t>> behind_indi;
		if (front_indi_num < pop.size()) {
			for (size_t i = 0; i < space_index.size(); ++i) {
				size_t temp = space_index[i];
				if (std::find(nondominate_space.begin(), nondominate_space.end(), temp) == nondominate_space.end()) {
					if (dominated_space.empty())
						dominated_space.emplace_back(temp);
					else if (std::find(dominated_space.begin(), dominated_space.end(), temp) == dominated_space.end()) {
						dominated_space.emplace_back(temp);
					}
				}
			}
			for (size_t i = 0; i < dominated_space.size(); ++i) {
				std::pair<size_t, std::vector<size_t>> temp;
				temp.first = dominated_space[i];
				for (size_t j = 0; j < space_index.size(); ++j) {
					if (space_index[j] == temp.first)
						temp.second.emplace_back(j);
				}
				behind_indi.insert(temp);
			}
		}
		space_all_indi.emplace_back(behind_indi);
		return space_all_indi;
	}
	
	pop_attach SPMOEA9::popAttach(const Population<Solution<>>& pop) {
		pop_attach all_indi;
		all_indi.emplace_back(spaceAttach(pop, m_mo_hlc->getObjspaceTree(), true));
		all_indi.emplace_back(spaceAttach(pop, m_mo_hlc->subspaceTree(), false));
		return all_indi;
	}

	std::vector<std::vector<size_t>> SPMOEA9::clusterExploitSpace(const KDTree& tree, const std::map<size_t,std::vector<size_t>>& frontspace) {
		//先取出子空间
		std::vector<size_t> front_space_inx;
		for (auto ss : frontspace) {
			front_space_inx.push_back(ss.first);
		}
		//再进行最优子空间的聚类
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> first_cluster;
		first_cluster.push_back(front_space_inx[0]);
		clustered.push_back(first_cluster);
		std::vector<size_t> select_flag(front_space_inx.size(),0);
		select_flag[0] = 1;
		size_t i = 0;
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			for (size_t j = 0; j < clustered[i].size(); ++j) {
				size_t inx = clustered[i][j];
				std::list<size_t> neighbors;
				tree.findNeighbor(inx, neighbors);
				for (size_t k = 0; k < front_space_inx.size(); ++k) {
					if (select_flag[k] == 0 && std::find(neighbors.begin(), neighbors.end(), front_space_inx[k]) != neighbors.end()) {
						clustered[i].push_back(front_space_inx[k]);
						select_flag[k] = 1;
					}
				}
			}
			if (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
				for (size_t k = 0; k < front_space_inx.size(); ++k) {
					if (select_flag[k] == 0) {
						std::vector<size_t> add_cluster;
						add_cluster.push_back(front_space_inx[k]);
						select_flag[k] = 1;
						clustered.push_back(add_cluster);
						i++;
						break;
					}
				}
			}
		}
		if (m_add_neighbor) {
			//再加入各类的一级邻域
			std::vector<size_t> flag(tree.size(), 0);//标记已选子空间
			for (size_t i = 0; i < front_space_inx.size(); ++i) {
				flag[front_space_inx[i]] = 1;
			}
			for (size_t i = 0; i < clustered.size(); ++i) {
				std::vector<size_t> m_neighbors;
				for (size_t j = 0; j < clustered[i].size(); ++j) {
					size_t inx = clustered[i][j];
					std::list<size_t> neighbors;
					tree.findNeighbor(inx, neighbors);
					for (auto k : neighbors) {
						if (m_neighbors.empty() || std::find(m_neighbors.begin(), m_neighbors.end(), k) == m_neighbors.end()) {
							m_neighbors.push_back(k);
						}
					}
				}
				for (auto j : m_neighbors) {
					if (flag[j] == 0 && std::find(clustered[i].begin(), clustered[i].end(), j) == clustered[i].end()) {
						clustered[i].push_back(j);
						flag[j] = 1;
					}
				}
			}
		}
		return clustered;
	}

	std::vector<size_t> SPMOEA9::getExploitSpace() {
		std::vector<size_t> exploit_space;
		//先找出子空间排序值为0的子空间
		std::vector<size_t> front_space;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank == 0) {
				front_space.push_back(i);
			}
		}
		for (size_t i = 0; i < front_space.size(); ++i) {
			exploit_space.push_back(front_space[i]);
		}
		if (m_add_neighbor) {
			//再看其他子空间是否为这些子空间的邻域
			for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
				for (size_t j = 0; j < front_space.size(); ++j) {
					std::list<size_t> neighbors;
					m_mo_hlc->subspaceTree().findNeighbor(front_space[j], neighbors);
					if (std::find(neighbors.begin(), neighbors.end(), i) != neighbors.end()) {
						if (std::find(exploit_space.begin(), exploit_space.end(), i) == exploit_space.end()) {
							exploit_space.push_back(i);
							break;
						}
					}
				}
			}
		}
		return exploit_space;
	}

	std::vector<size_t> SPMOEA9::getExploreSpace() {
		//先找出子空间排序值为0的子空间
		std::vector<size_t> front_space;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank == 0) {
				front_space.push_back(i);
			}
		}
		//先将非前沿子空间加入
		std::vector<size_t> explore_space;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (std::find(front_space.begin(), front_space.end(), i) == front_space.end()) {
				explore_space.push_back(i);
			}
		}

		if (m_add_neighbor) {
			//将已加入的子空间中的前排空间的邻域去掉
			std::vector<size_t> select_flag(explore_space.size(),0);
			for (size_t i = 0; i < explore_space.size(); ++i) {
				for (size_t j = 0; j < front_space.size(); ++j) {
					std::list<size_t> neighbors;
					m_mo_hlc->subspaceTree().findNeighbor(front_space[j], neighbors);
					if (std::find(neighbors.begin(), neighbors.end(), explore_space[i]) == neighbors.end()) {
						if (j == front_space.size() - 1) {
							select_flag[i]=1;
						}
					}
					else {
						break;
					}
				}
			}
			std::vector<size_t> select_index;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 1) {
					select_index.push_back(explore_space[i]);
				}
			}
			explore_space.clear();
			explore_space = select_index;
		}
		return explore_space;
	}

	std::vector<std::vector<size_t>> SPMOEA9::clusterExploreSpace(const KDTree& tree) {
		//先得到探索子空间索引
		std::vector<size_t> explore_spaces = getExploreSpace();
		std::vector<std::vector<size_t>> explore_cluster;
		//再根据探索潜力值聚类
		//explore_cluster = getMO_HLC().clusterSubspace(explore_spaces,tree.size(), m_add_neighbor);
		//直接作为一类
		explore_cluster.push_back(explore_spaces);
		return explore_cluster;
	}

	std::vector<size_t> SPMOEA9::assignBasinResource(const std::vector<std::vector<size_t>>& potential_inds) {
		//归一化的潜力值
		std::vector<Real> resource;
		size_t total_num = 0;
		Real sum = 0.;
		/*for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
			auto basin = m_mo_hlc->getBasinInfo(i);
			if (basin.flag == "explore") {
				resource.push_back(basin.m_basin_explore_potential);
			}
			else {
				resource.push_back(basin.m_basin_exploit_potential);
			}
			sum += resource.back();
		}*/
		for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
			auto basin = m_mo_hlc->getBasinInfo(i);
			resource.push_back(basin.m_basin_potential);
			sum += resource.back();
		}
		for (size_t i = 0; i < potential_inds.size(); ++i) {
			Real temp_potential = 0.1;
			resource.push_back(temp_potential);
			sum += temp_potential;
		}
		for (auto& i : resource) {
			i = i / sum;
		}
		std::vector<size_t> assign_resource;
		for (size_t i = 0; i < resource.size(); ++i) {
			auto tmp = std::ceil(m_pop->size() * resource[i]);
			assign_resource.push_back(tmp);
			total_num += tmp;
		}
		//去除多余个体
		while (total_num > m_pop->size()) {
			size_t max = 0;
			size_t inx=0;
			for (size_t i = 0; i < resource.size(); ++i) {
				if (assign_resource[i] > max) {
					max = assign_resource[i];
					inx = i;
				}
			}
			assign_resource[inx] = assign_resource[inx] - 1;
			total_num -= 1;
		}
		return assign_resource;
	}

	void SPMOEA9::recordMetrics(Problem *pro) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		if (m_pop->iteration() < 1) {
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < m_pop->size(); ++i)
				objs.emplace_back(&m_pop->at(i).objective());
			std::vector<int> rank;
			ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_pop->at(i).setFitness(rank[i]);
				if (rank[i] == 0)
					m_history_front_sols.emplace_back(std::make_shared<Solution<>>(m_pop->at(i)));
			}
			//updateObjSpaceInfo(m_pop->getParent(), pro);//更新子空间信息
		}
		std::vector<Solution<>> temp_pop;
		size_t count1,count2;
		count1 = count2 = 0;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			if (m_pop->at(i).fitness() == 0) {
				temp_pop.emplace_back(m_pop->at(i));
				count1++;
			}
			else if (m_pop->at(i).fitness() == 1) {
				count2++;
			}
		}
		m_IGD = pro->optimaBase()->invertGenDist(temp_pop);
		m_R1 = (Real)count1 / m_pop->size();
		m_R2 = (Real)count1 / (1 + count2)/ m_pop->size();
		m_R3 = (Real)count1 / (1 + m_pop->size()-count1-count2)/ m_pop->size();
		record();//store metrics data
	}

	void SPMOEA9::recordHisFront(Problem *pro) {
		//记录历史前沿
		/**********************/
		/*   将前沿个体写入文本
		/**********************/
		size_t num_eval = evaluations();
		/*if (num_eval == 2600) {
			size_t a = 1;
		}*/
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		if ((num_eval == 300000 && num_var < 5) || (num_eval == 500000 && num_var >= 5 && num_var<10) || (num_eval == 700000 && num_var >=10)) {
			std::string pro_name = CAST_CONOP(pro)->name();
			std::string file_name2 = "E:/Gitlab/OFEC/result/optima_samples/" + std::to_string(num_obj) + "_objs/" \
				+ std::to_string(num_var) + "_vars/" + pro_name + ".txt";
			std::ofstream out2(file_name2);
			if (out2) {
				for (size_t i = 0; i < num_obj + num_var; ++i) {
					out2 << num_obj << " ";
				}
				out2 << std::endl;
				for (size_t i = 0; i < num_obj + num_var; ++i) {
					out2 << num_var << " ";
				}
				out2 << std::endl;
				for (auto& me : m_history_front_sols) {
					auto var = me->variable().vect();
					auto obj = me->objective();
					for (size_t i = 0; i < var.size(); ++i) {
						out2 << var[i] << " ";
					}
					for (size_t i = 0; i < obj.size(); ++i) {
						out2 << obj[i] << " ";
					}
					out2 << std::endl;
				}
				out2.close();
			}
		}
	}

	void SPMOEA9::updatePopAge() {
		for (size_t i = 0; i < m_pop->size(); i++) {
			size_t temp = m_pop->at(i).surviveAge();
			m_pop->at(i).setCounter(temp + 1);
		}
	}

	int SPMOEA9_pop::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		return 0;
	}
}
