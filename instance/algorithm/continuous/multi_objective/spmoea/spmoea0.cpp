#include "spmoea0.h"
#include<algorithm>
#include"../../../../../utility/clustering/dbscan.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA0::initialize_() {
		SPMOEA::initialize_();
		auto& v = *m_param;
		//m_add_neighbor = v.get<bool>(("cluster add neighbors"));
		//m_sample_in_basin = v.get<bool>(("sample in basin"));
		//m_subspace_select = v.get<bool>(("select by subspace"));
		//m_evolve_by_potential = v.get<bool>(("evolve by potential"));
		//m_mean_potential = v.get<bool>(("mean basin potential"));
		m_exploit_pop_size=v.get<int>("population size");
		m_explore_pop_size = v.get<int>("population size");
		m_num_e_e.resize(2);
	}

	void SPMOEA0::run_() {
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

	void SPMOEA0::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = CAST_CONOP(m_problem.get())->getOptima().invertGenDist(*m_pop);
		entry.push_back(m_IGD);
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA0::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[0].push_back(&getPop()[i][j]);
				}
				/*for (size_t i = 0; i < m_pop->getOffspring().size(); ++i)
					m_solution[1].push_back(&m_pop->getOffspring()[i]);*/
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA0::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
        //先在子空间预采样，根据比较结果生成多种群
		size_t num = getMO_HLC().numSubspace();
		SPMOEA_pop temp_pop(num, pro);
		for (size_t i = 0; i < num; ++i) {
			auto bound = getMO_HLC().subspaceTree().getBox(i);
			for (size_t k = 0; k < bound.size(); ++k) {
				temp_pop[i].variable().vect()[k] = bound[k].first + 1. / 2 * (bound[k].second - bound[k].first);
			}
			temp_pop[i].evaluate(pro,alg);
		}
		NDSort(temp_pop);
		//更新子空间信息
		updateVarSpaceInfo(temp_pop, pro, rnd);
		updateHistoryInfo(temp_pop, pro);
		auto space_clusters = clusterSubspaces();
		auto front_ssp = getFrontSubspace();
		setPreFrontSubspace(front_ssp);
		generatePop(space_clusters,pro,alg,rnd);


		SPMOEA::updateArchive(archiveNum(),pro);
		SPMOEA::updateObjRange(temp_pop,pro);
		updateObjSpace();
		//还原所有子空间采样状态
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).add_flag = "no";
		}
	}

	void SPMOEA0::generatePop(const std::vector<std::vector<size_t>>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd) {
		for (size_t i = 0; i < space_clusters.size(); ++i) {
			size_t num = space_clusters[i].size();
			size_t pop_size;
            if(i == space_clusters.size() - 1) {
				pop_size = m_explore_pop_size;
			}
			else {
				pop_size = m_exploit_pop_size;
			}
			SPMOEA_pop temp_pop(pop_size, pro);
			size_t count = 0;
			
			while (count < pop_size) {
				size_t inx = std::floor(num * rnd->uniform.next());
				auto bound = getMO_HLC().subspaceTree().getBox(space_clusters[i][inx]);
				if (i == space_clusters.size() - 1) {
					//探索子空间的采样方式:均匀，邻域采样
					bool if_random = true;
					if (if_random) {
						for (size_t k = 0; k < bound.size(); ++k) {
							temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
						}
					}
					else {
						//在该子空间内的前沿个体附近进行局部采样
						auto& present_ind = getMO_HLC().getSubspaceInfo(space_clusters[i][inx]).m_subspace_front_sol;
						if (present_ind.size() == 0) {
							for (size_t k = 0; k < bound.size(); ++k) {
								temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
							}
						}
						else if (present_ind.size() == 1) {
							auto pos = present_ind[0]->variable().vect();
							auto& his_ind = getMO_HLC().getSubspaceInfo(space_clusters[i][inx]).m_history_inds;
							if (his_ind.size() > 1) {
								//在历史解中随机选一个解
								size_t ind_inx;
								std::vector<Real> start_pos;
								do {
									ind_inx = std::floor(his_ind.size() * rnd->uniform.next());
									start_pos = his_ind[ind_inx]->variable().vect();
								} while (!ifSame(pos, start_pos));
								std::vector<Real> new_pos(bound.size(),0.);
								Real ratio =rnd->uniform.next();
								for (size_t k = 0; k < bound.size(); ++k) {
									new_pos[k] = (1+ratio) * (pos[k]-start_pos[k]) + start_pos[k];
								}
								repairSol(new_pos, bound, rnd);
								temp_pop[count].variable().vect() = new_pos;
							}
							else {
								for (size_t k = 0; k < bound.size(); ++k) {
									temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
								}
							}
						}
						else {
							PopSBX<> temp_pop2(present_ind.size(), pro);
							for (size_t k = 0; k < present_ind.size(); ++k) {
								temp_pop2[k].variable() = present_ind[k]->variable();
								temp_pop2[k].objective() = present_ind[k]->objective();
							}
							Solution<> temp_sol1(temp_pop2[0]);
							Solution<> temp_sol2(temp_pop2[1]);
							std::vector<size_t> p(2);
							do {//不重复采样
								do {
									p[0] = temp_pop2.tournamentSelection(pro, rnd);
									p[1] = temp_pop2.tournamentSelection(pro, rnd);
								} while (p[1] == p[0]);
								temp_pop2.crossover(p[0], p[1], temp_sol1, temp_sol2, pro, rnd);
								temp_pop2.mutate(temp_sol1, pro, rnd);
								//子代越界处理,允许探索到邻域？
								repairSol(temp_sol1.variable().vect(), bound, rnd);
							} while (temp_sol1.variable().vect() == temp_pop2[p[0]].variable().vect() || \
								temp_sol1.variable().vect() == temp_pop2[p[1]].variable().vect());
							
							temp_pop[count].variable() = temp_sol1.variable();
						}
					}
				}
				else {
					for (size_t k = 0; k < bound.size(); ++k) {
						temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
					}
				}
				count++;
				getMO_HLC().getSubspaceInfo(space_clusters[i][inx]).add_flag = "yes";
			}
			temp_pop.evaluate(pro, alg);
			NDSort(temp_pop);
			//初始化子代
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring()[j] = temp_pop[j];
				temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			}
			/*temp_pop.setRate(getCr(), getMr());
			temp_pop.setEta(getCeta(), getMeta());*/
			if (i == space_clusters.size() - 1) {
				temp_pop.setPopState("explore");
				for (size_t j = 0; j < space_clusters[i].size(); ++j) {
					getMO_HLC().getSubspaceInfo(space_clusters[i][j]).search_flag = "explore";
				}
			}
			else {
				temp_pop.setPopState("exploit");
				for (size_t j = 0; j < space_clusters[i].size(); ++j) {
					getMO_HLC().getSubspaceInfo(space_clusters[i][j]).search_flag = "exploit";
				}
			}
			temp_pop.setLocateSubspace(space_clusters[i]);
			//更新子空间信息
			updateVarSpaceInfo(temp_pop, pro, rnd);
			//更新历史前沿
			updateHistoryInfo(temp_pop, pro);
			getPop().append(temp_pop);

			////输出种群
			//for (size_t j = 0; j < temp_pop.size(); ++j) {
			//	std::cout << temp_pop[j].variable()[0] << "      " << temp_pop[j].variable()[1] << std::endl;
			//}
			//std::cout << "*******************" << std::endl;
		}
	}

	void SPMOEA0::generateGoodPop(const std::vector<size_t>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd) {
		//先得到子空间当前最优前沿解，
		std::vector<std::shared_ptr<Solution<>>> pop;
		for (size_t i = 0; i < space_clusters.size(); ++i) {
			auto space_inx = space_clusters[i];
			for (size_t j = 0; j < getMO_HLC().getSubspaceInfo(space_inx).m_represent_sol.size(); ++j) {
				pop.emplace_back(getMO_HLC().getSubspaceInfo(space_inx).m_represent_sol[j]);
			}
			getMO_HLC().getSubspaceInfo(space_clusters[i]).search_flag = "exploit";
		}
		//再根据最优前沿解生成子种群
		size_t pop_size=m_exploit_pop_size;;
		SPMOEA_pop temp_pop(pop_size, pro);
		if (pop_size < pop.size()) {
			std::vector<size_t> add_inx;
			for (size_t i = 0; i < pop_size; ++i) {
				auto inx = std::floor(pop.size()*rnd->uniform.next());
				if (add_inx.empty() || std::find(add_inx.begin(), add_inx.end(), inx) == add_inx.end()) {
					temp_pop[i] = *pop[inx];
				}
			}
		}
		else {
			for (size_t i = 0; i < pop.size(); ++i) {
				temp_pop[i] = *pop[i];
			}
			//额外加入部分个体
			SPMOEA_pop temp_pop2(pop_size-pop.size(), pro);
			size_t count = 0;
			while (count < pop_size - pop.size()) {
				size_t inx = std::floor(space_clusters.size() * rnd->uniform.next());
				auto bound = getMO_HLC().subspaceTree().getBox(space_clusters[inx]);
				for (size_t k = 0; k < bound.size(); ++k) {
					temp_pop2[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
				}
				count++;
				getMO_HLC().getSubspaceInfo(space_clusters[inx]).add_flag = "yes";
			}
			temp_pop2.evaluate(pro, alg);
			for (size_t i = 0; i < temp_pop2.size(); ++i) {
				temp_pop[pop.size() + i] = temp_pop2[i];
			}
			NDSort(temp_pop2);
			//更新子空间信息
			updateVarSpaceInfo(temp_pop2, pro, rnd);
			//更新历史前沿
			updateHistoryInfo(temp_pop2, pro);
		}
		//初始化子代
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
		}
		/*temp_pop.setRate(getCr(), getMr());
		temp_pop.setEta(getCeta(), getMeta());*/
		temp_pop.setPopState("exploit");
		temp_pop.setLocateSubspace(space_clusters);

		getPop().append(temp_pop);
	}

	void SPMOEA0::generateOffspring(Problem *pro, Random *rnd) {
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	auto spaces = getPop()[i].getLocateSubspace();
		//	if (getPop()[i].getPopState() == "exploit") {
		//		for (size_t j = 0; j < getPop()[i].size(); j += 2) {
		//			std::vector<size_t> p(2);
		//			do {//不重复采样
		//				do {
		//					p[0] = getPop()[i].tournamentSelection(pro, rnd);
		//					p[1] = getPop()[i].tournamentSelection(pro, rnd);
		//				} while (p[1] == p[0]);
		//				getPop()[i].crossover(p[0], p[1], getPop()[i].getOffspring()[j], getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				//
		//				std::vector<std::vector<Real>> ppp;
		//				ppp.emplace_back(getPop()[i][p[0]].variable().vect());
		//				ppp.emplace_back(getPop()[i][p[1]].variable().vect());
		//				ppp.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
		//				ppp.emplace_back(getPop()[i].getOffspring()[j+1].variable().vect());
		//				// 
		//				//子代越界处理,允许探索到邻域？
		//				repairSol(getPop()[i].getOffspring()[j].variable().vect(), spaces, rnd);
		//				repairSol(getPop()[i].getOffspring()[j + 1].variable().vect(), spaces, rnd);
		//			} while (getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[1]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[1]].variable().vect());

		//			getPop()[i].getOffspring()[j].setTimeEvaluate(0);
		//			getPop()[i].getOffspring()[j + 1].setTimeEvaluate(0);
		//			for (size_t k = 0; k < spaces.size(); ++k) {
		//				getMO_HLC().getSubspaceInfo(spaces[k]).add_flag = "yes";
		//			}
		//		}
		//	}
		//	else {
		//		//随机采样，锦标赛采样，以及基于覆盖率的采样相结合
		//		//根据探索子空间的覆盖率确定采样子空间
		//		//探索种群在所属的子空间中随机采样还是在已知代表解附近采样
		//		bool if_random = true;
		//		if (if_random) {
		//			for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//				size_t inx = (size_t)std::floor(rnd->uniform.next() * spaces.size());
		//				auto bound = getMO_HLC().subspaceTree().getBox(spaces[inx]);
		//				std::vector<Real> temp_sol;
		//				for (size_t k = 0; k < bound.size(); ++k) {
		//					temp_sol.push_back(bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first));
		//				}
		//				getPop()[i].getOffspring()[j].variable().vect() = temp_sol;
		//				getPop()[i].getOffspring()[j].setTimeEvaluate(0);
		//				getMO_HLC().getSubspaceInfo(spaces[inx]).add_flag = "yes";
		//			}
		//		}
		//		else {
		//			//在该子空间内的前沿个体附近进行局部采样
		//			for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//				size_t inx = (size_t)std::floor(rnd->uniform.next() * spaces.size());
		//				auto bound = getMO_HLC().subspaceTree().getBox(spaces[inx]);

		//				auto& present_ind = getMO_HLC().getSubspaceInfo(spaces[inx]).m_subspace_front_sol;
		//				std::vector<Real> temp_sol(bound.size(),0.);
		//				if (present_ind.size() == 0) {
		//					for (size_t k = 0; k < bound.size(); ++k) {
		//						temp_sol[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
		//					}
		//				}
		//				else if (present_ind.size() == 1) {
		//					auto pos = present_ind[0]->variable().vect();
		//					auto& his_ind = getMO_HLC().getSubspaceInfo(spaces[inx]).m_history_inds;
		//					if (his_ind.size() > 1) {
		//						//在历史解中随机选一个解
		//						size_t ind_inx;
		//						std::vector<Real> start_pos;
		//						do {
		//							ind_inx = std::floor(his_ind.size() * rnd->uniform.next());
		//							start_pos = his_ind[ind_inx]->variable().vect();
		//						} while (!ifSame(pos, start_pos));
		//						std::vector<Real> new_pos(bound.size(), 0.);
		//						Real ratio = rnd->uniform.next();
		//						for (size_t k = 0; k < bound.size(); ++k) {
		//							new_pos[k] = (1 + ratio) * (pos[k] - start_pos[k]) + start_pos[k];
		//						}
		//						repairSol(new_pos, bound, rnd);
		//						temp_sol = new_pos;
		//					}
		//					else {
		//						for (size_t k = 0; k < bound.size(); ++k) {
		//							temp_sol[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
		//						}
		//					}
		//				}
		//				else {
		//					SPMOEA_pop temp_pop2(present_ind.size(), pro);
		//					for (size_t k = 0; k < present_ind.size(); ++k) {
		//						temp_pop2[k].variable() = present_ind[k]->variable();
		//						temp_pop2[k].objective() = present_ind[k]->objective();
		//					}
		//					std::vector<size_t> p(2);
		//					do {//不重复采样
		//						do {
		//							p[0] = temp_pop2.tournamentSelection(pro, rnd);
		//							p[1] = temp_pop2.tournamentSelection(pro, rnd);
		//						} while (p[1] == p[0]);
		//						temp_pop2.crossover(p[0], p[1], temp_pop2.getOffspring()[0], temp_pop2.getOffspring()[1], pro, rnd);
		//						temp_pop2.mutate(temp_pop2.getOffspring()[0], pro, rnd);
		//						//
		//						std::vector<std::vector<Real>> ppp;
		//						ppp.emplace_back(temp_pop2[p[0]].variable().vect());
		//						ppp.emplace_back(temp_pop2[p[1]].variable().vect());
		//						ppp.emplace_back(temp_pop2.getOffspring()[0].variable().vect());
		//						//
		//						//子代越界处理,允许探索到邻域？
		//						repairSol(temp_pop2.getOffspring()[0].variable().vect(), bound, rnd);
		//					} while (temp_pop2.getOffspring()[0].variable().vect() == temp_pop2[p[0]].variable().vect() || \
		//						temp_pop2.getOffspring()[0].variable().vect() == temp_pop2[p[1]].variable().vect());

		//					temp_sol = temp_pop2.getOffspring()[0].variable().vect();
		//				}
		//				getPop()[i].getOffspring()[j].variable().vect() = temp_sol;
		//				getPop()[i].getOffspring()[j].setCounter(0);
		//				getMO_HLC().getSubspaceInfo(spaces[inx]).add_flag = "yes";
		//			}
		//		}
		//	}
		//}
	}

	int SPMOEA0::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		//auto& v = *m_param;
		//m_evolve_by_potential = v.get<bool>("evolve by potential");
		/*
		根据子空间的explore还是exploit状态，决定种群采样方式，通过局部前沿检测方式动态调整子空间状态
		*/
		/******************************************
		     对开发子空间，更新前沿子空间聚类，
		******************************************/
		/*************************************************************************
		 重新聚类子空间，对于之前在前沿的，继续搜索，对于新加入的，初始化一个种群
		**************************************************************************/
		auto space_clusters = clusterSubspaces();
		auto pre_front = getPreFrontSubspace();
		auto front_subspace = getFrontSubspace();
		auto behind_subspace = space_clusters.back();

		if (!ifSame(pre_front, front_subspace)) {
			updateMultiPop(space_clusters, pre_front, pro, alg, rnd);
			setPreFrontSubspace(front_subspace);
		}
		/*****************************************************
			 对探索子空间，发现子空间中是否具有局部前沿，
		*****************************************************/
		std::vector<size_t> local_pareto_subspaces;//具有局部前沿的子空间索引
		std::vector<size_t> degenerate_subspaces;//从具有局部前沿退化到没有局部前沿的子空间索引
		for (size_t i = 0; i < behind_subspace.size(); ++i) {
			//根据子空间的采用数据分布，判断是否具有局部前沿
			if (getMO_HLC().getSubspaceInfo(behind_subspace[i]).add_flag == "yes") {//子空间有新数据
				if (isLocalParetoRegion(behind_subspace[i],pro,alg,rnd)) {
					local_pareto_subspaces.push_back(behind_subspace[i]);
				}
				else if(getMO_HLC().getSubspaceInfo(behind_subspace[i]).search_flag=="exploit") {
					degenerate_subspaces.push_back(behind_subspace[i]);
				}
			}
		}
		//还原所有子空间采样状态
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).add_flag = "no";
		}
		if (!degenerate_subspaces.empty()) {
			//删除该子种群
			for (size_t i = 0; i < degenerate_subspaces.size(); ++i) {
				size_t pop_inx = 0;
				for (size_t j = 0; j < getPop().size(); ++j) {
					if (getPop()[j].getLocateSubspace()[0] == degenerate_subspaces[i]) {
						pop_inx = j;
						break;
					}
				}
				getPop().remove(getPop().begin() + pop_inx);
			}
		}
		if (!degenerate_subspaces.empty() || !local_pareto_subspaces.empty()) {
			//先删除探索子种群
			bool has_explore_pop = false;
			for (size_t i = 0; i < getPop().size(); ++i) {
				if (getPop()[i].getPopState() == "explore") {
					has_explore_pop = true;
					break;
				}
			}
			if (has_explore_pop) {
				getPop().remove(getPop().end()-1);
			}
		}
		if (!local_pareto_subspaces.empty()) {
			//对这些子空间聚类，并在每个类中生成子种群搜索，并且将聚类的子空间加入总的类中，并去除对应的探索子空间
			//除了邻域关系之外，还要看两个子空间的前沿能否聚类
			//如果判断是局部子空间，还要判断上一代的种群状态是什么？或者子空间的状态是什么？探索还是开发
			for (size_t i = 0; i < local_pareto_subspaces.size(); ++i) {
				std::vector<size_t> temp_space;
				temp_space.push_back(local_pareto_subspaces[i]);
				//先看现有子种群有没有在该空间演化，再新增
				bool pop_exist = false;
				for (size_t j = 0; j < getPop().size(); ++j) {
					auto temp = getPop()[j].getLocateSubspace();
					if (std::find(temp.begin(), temp.end(), local_pareto_subspaces[i]) != temp.end()) {
						pop_exist = true;
						break;
					}
				}
				if (!pop_exist) {
					generateGoodPop(temp_space, pro, alg, rnd);
				}
				//从后排子空间中删除该子空间
				size_t inx = std::find(behind_subspace.begin(),behind_subspace.end(),local_pareto_subspaces[i])- behind_subspace.begin();
				behind_subspace.erase(behind_subspace.begin()+inx);
			}
		}
		bool has_explore_pop = false;
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "explore") {
				has_explore_pop = true;
				break;
			}
		}
		if (!has_explore_pop) {
			//在后排子空间生成探索子种群
			if (!behind_subspace.empty()) {
				std::vector<std::vector<size_t>> temp_space;
				temp_space.emplace_back(behind_subspace);
				generatePop(temp_space, pro, alg, rnd);
			}
		}
		/*****************************************************
			      产生子代,采样，评价，更新子空间信息
				  子代采样的目的：产生新的好解
		*****************************************************/
		generateOffspring(pro, rnd);//子代产生容许采样到邻域
		//评价子代，更新子空间信息
		int tag = EvaluationTag::kNormalEval;
		size_t ind_counter = 0;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); j++) {
				//评价子代
				tag = getPop()[i].getOffspring()[j].evaluate(pro, alg);
				if (tag != EvaluationTag::kNormalEval)
					break;
				ind_counter++;
			}
			if (tag != EvaluationTag::kNormalEval)
				break;
		}
		if (tag != EvaluationTag::kNormalEval) {
			return tag;
		}
		SPMOEA_pop offspring_pop(ind_counter,pro);
		size_t count = 0;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); j++) {
				offspring_pop[count + j] = getPop()[i].getOffspring()[j];
			}
			count += getPop()[i].size();
		}
		//更新子空间信息
		NDSort(offspring_pop);
		updateVarSpaceInfo(offspring_pop, pro, rnd);
		//更新历史信息
		updateHistoryInfo(offspring_pop, pro);

		/*****************************************************
				           更新子种群
				环境选择的目的：更有利于产生新的好解
		*****************************************************/
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			getPop()[i].envirSelection(getPop()[i], getPop()[i].getOffspring(),pro);
			
			//if (i == getPop().size() - 1) {
			//	//输出种群
			//	for (size_t j = 0; j < getPop()[i].size(); ++j) {
			//		std::cout << getPop()[i][j].variable()[0] << "      " << getPop()[i][j].variable()[1] << std::endl;
			//	}
			//	std::cout << "************************" << std::endl;
			//}
			
			//getPop()[i].iteration();
		}

		SPMOEA::updateObjRange(offspring_pop,pro);

		//更新存档
		SPMOEA::updateArchive(archiveNum(), pro);
		//updateObjSpace();

		SPMOEA::recordMetrics(pro, alg);

		return tag;
	}

	bool SPMOEA0::isLocalParetoRegion(size_t inx,Problem *pro,Algorithm *alg,Random *rnd) {
		//根据子空间的前沿解，找出两个分解超平面
		//找一个支配它的最近的2个点，并生成一个对称点，评价对称点，看是否被子空间前沿支配
		//或者找出前沿个体张成空间的角点，比较角点与前沿解的支配关系
		//先得到子空间的前沿个体
		bool islocal = false;
		
		// method 1
		auto &present_ind = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		if (present_ind.size() < num_var) {
			islocal = false;
		}
		else {
			std::vector<Real> dim_span;
			std::map<size_t, std::pair<Real, Real>> dim_slice;
			std::vector<size_t> lower_inx;
			std::vector<size_t> upper_inx;
			auto& his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
			//寻找每一维的切面
			for (size_t i = 0; i < num_var; ++i) {
				Real min_v = (Real)INT16_MAX;
				Real max_v = (Real)INT16_MAX * -1;
				for (size_t j = 0; j < present_ind.size(); ++j) {
					if (min_v > present_ind[j]->variable().vect()[i]) {
						min_v = present_ind[j]->variable().vect()[i];
					}
					if (max_v < present_ind[j]->variable().vect()[i]) {
						max_v = present_ind[j]->variable().vect()[i];
					}
				}
				dim_span.push_back(max_v - min_v);
				dim_slice.insert(std::make_pair(i, std::make_pair(min_v, max_v)));
			}
			//找出最小的跨度的维
			auto min_value = std::min_element(dim_span.begin(), dim_span.end());
			auto min_dim_inx = std::distance(dim_span.begin(), min_value);
			//封闭超矩形的中心坐标
			std::vector<Real> center;
			for (size_t i = 0; i < num_var; ++i) {
				center.push_back((dim_slice[i].first + dim_slice[i].second) / 2);
			}
			//将中心点沿着最小跨度维平移得到2个测试点
			Real shift_v = 0.;
			Real d1 = dim_slice[min_dim_inx].second - center[min_dim_inx];
			//Real d2 = center[min_dim_inx]- dim_slice[min_dim_inx].first;
			Real d3 = getMO_HLC().subspaceTree().getBox(inx)[min_dim_inx].second - center[min_dim_inx];
			Real d4 = center[min_dim_inx] - getMO_HLC().subspaceTree().getBox(inx)[min_dim_inx].first;
			Real max_b = d3 < d4 ? d3 : d4;
			shift_v = d1 + rnd->uniform.next() * (max_b - d1);
			std::vector<Real> ind1, ind2;
			ind1 = ind2 = center;
			ind1[min_dim_inx] += shift_v;
			ind2[min_dim_inx] -= shift_v;
			SPMOEA_pop temp_pop(2, pro);
			temp_pop[0].variable().vect() = ind1;
			temp_pop[1].variable().vect() = ind2;
			temp_pop.evaluate(pro, alg);
			getMO_HLC().getSubspaceInfo(inx).add_flag = "yes";
			NDSort(temp_pop);
			//判断两个解是否被当前子空间前沿支配
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				bool flag = false;
				for (size_t j = 0; j < present_ind.size(); ++j) {
					auto relation = temp_pop[i].compare(*present_ind[j], CAST_CONOP(pro)->optimizeMode());
					if (relation == Dominance::kDominated) {
						flag = true;
					}
				}
				if (!flag) {
					break;
				}
				if (i == temp_pop.size() - 1) {
					islocal = true;
				}
			}
			//更新子空间信息
			updateVarSpaceInfo(temp_pop, pro, rnd);
			//更新历史前沿
			updateHistoryInfo(temp_pop, pro);
		}
		return islocal;
		
        ////method 2
		////找出最小的跨度的维
		//auto &present_ind = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		//size_t num_var = CAST_CONOP(pro)->numberVariables();
		//if (present_ind.size() < num_var) {
		//	islocal = false;
		//}
		//else {
		//	std::vector<Real> dim_span;
		//	std::map<size_t, std::pair<Real, Real>> dim_slice;
		//	std::vector<size_t> lower_inx;
		//	std::vector<size_t> upper_inx;
		//	auto& his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		//	//寻找每一维的切面
		//	for (size_t i = 0; i < num_var; ++i) {
		//		Real min_v = (Real)INT16_MAX;
		//		Real max_v = (Real)INT16_MAX * -1;
		//		for (size_t j = 0; j < present_ind.size(); ++j) {
		//			if (min_v > present_ind[j]->variable().vect()[i]) {
		//				min_v = present_ind[j]->variable().vect()[i];
		//			}
		//			if (max_v < present_ind[j]->variable().vect()[i]) {
		//				max_v = present_ind[j]->variable().vect()[i];
		//			}
		//		}
		//		dim_span.push_back(max_v - min_v);
		//		dim_slice.insert(std::make_pair(i, std::make_pair(min_v, max_v)));
		//	}
		//	//找出最小的跨度的维
		//	auto min_value = std::min_element(dim_span.begin(), dim_span.end());
		//	auto min_dim_inx = std::distance(dim_span.begin(), min_value);

		//	for (size_t i = 0; i < his_ind.size(); ++i) {
		//		auto temp_pos = his_ind[i]->variable().vect();
		//		if (temp_pos[min_dim_inx] < dim_slice[min_dim_inx].first) {
		//			lower_inx.push_back(i);
		//		}
		//		if (temp_pos[min_dim_inx] > dim_slice[min_dim_inx].second) {
		//			upper_inx.push_back(i);
		//		}
		//	}
		//	if (lower_inx.empty() || upper_inx.empty()) {
		//		islocal=false;
		//	}
		//	else {
		//		islocal=true;
		//	}
		//}
		//return islocal;

  //      //method 3
		//auto &present_ind = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		//size_t num_var = CAST_CONOP(pro)->numberVariables();
		//if (present_ind.size() < num_var) {
		//	islocal = false;
		//}
		//else {
		//	islocal = true;
		//}
		//if (islocal) {
		//	/*********************************************************
		//	*                   根据样本构造训练数据 
		//	*********************************************************/
		//	std::vector<std::vector<Real>> train_data;
		//	std::vector<std::vector<Real>> origin_data;
		//	for (size_t i = 0; i < present_ind.size(); ++i) {
		//		auto temp_pos = present_ind[i]->variable().vect();
		//		origin_data.emplace_back(temp_pos);
		//	}
		//	size_t sample_num = 10 * num_var;
		//	bool transform_coordinate = true;//是否进行坐标变换
		//	if (present_ind.size() <= sample_num) {
		//		for (size_t i = 0; i < present_ind.size(); ++i) {
		//			auto temp_pos = present_ind[i]->variable().vect();
		//			train_data.push_back(temp_pos);
		//		}
		//	}
		//	else {
		//		std::vector<size_t> rand_seq;
		//		rand_seq.resize(present_ind.size());
		//		for (size_t i = 0; i < present_ind.size(); ++i) {
		//			rand_seq[i] = i;
		//		}
		//		rnd->uniform.shuffle(rand_seq.begin(), rand_seq.end());
		//		for (size_t i = 0; i < sample_num; ++i) {
		//			auto temp_pos = present_ind[rand_seq[i]]->variable().vect();
		//			train_data.push_back(temp_pos);
		//		}
		//	}
		//	std::vector<Real> t_origin(num_var, 0.);
		//	if (transform_coordinate) {
		//		//离哪一维的平均距离最近，然后其余维取大值，得到参考坐标原点
		//		auto bound = getMO_HLC().subspaceTree().getBox(inx);
		//		std::vector<std::pair<Real, Real>> dim_dist;
		//		for (size_t j = 0; j < num_var; ++j) {
		//			Real sum = 0.;
		//			for (size_t k = 0; k < train_data.size(); ++k) {
		//				sum += train_data[k][j];
		//			}
		//			sum /= train_data.size();
		//			dim_dist.emplace_back(std::make_pair(std::fabs(sum - bound[j].first), std::fabs(sum - bound[j].second)));
		//		}
		//		std::vector<Real> all_dist;
		//		for (size_t i = 0; i < dim_dist.size(); ++i) {
		//			all_dist.push_back(dim_dist[i].first);
		//			all_dist.push_back(dim_dist[i].second);
		//		}
		//		auto min_value = std::min_element(all_dist.begin(), all_dist.end());
		//		auto min_dim_inx = std::distance(all_dist.begin(), min_value);
		//		size_t select_dim = min_dim_inx / 2;
		//		//构建原点坐标
		//		for (size_t i = 0; i < num_var; ++i) {
		//		    t_origin[i] = dim_dist[i].first < dim_dist[i].second ? bound[i].first : bound[i].second;
		//		}
		//		//然后根据参考坐标原点变换训练数据和原始数据
		//		for (size_t i = 0; i < train_data.size(); ++i) {
		//			for (size_t j = 0; j < num_var; ++j) {
		//				train_data[i][j] -= t_origin[j];
		//			}
		//		}
		//		for (size_t i = 0; i < origin_data.size(); ++i) {
		//			for (size_t j = 0; j < num_var; ++j) {
		//				origin_data[i][j] -= t_origin[j];
		//			}
		//		}
		//	}
		//	for (size_t i = 0; i < train_data.size(); ++i) {
		//		//在倒数第二个位置插入一个常数值
		//		train_data[i].insert(train_data[i].end() - 1, 1.);
		//	}

		//	/*********************************************************
		//	*                       训练拟合平面
		//	*********************************************************/
		//	Real learn_ratio =1;
		//	size_t iterate_times = 10000;
		//	auto weigh = linearRegression(train_data, learn_ratio, iterate_times);
		//	//平移超平面直到覆盖子空间前沿个体
		//	auto upper_weigh = weigh;
		//	auto lower_weigh = weigh;
		//	auto offset1=shiftPlane(origin_data, upper_weigh, 1);
		//	upper_weigh.back() -= offset1;
		//	auto offset2 = shiftPlane(origin_data, lower_weigh, 0);
		//	lower_weigh.back() -= offset2;
		//	
		//	//找出上下超平面外的个体
		//	auto& his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		//	std::vector<size_t> lower_inx;
		//	std::vector<size_t> upper_inx;
		//	for (size_t i = 0; i < his_ind.size(); ++i) {
		//		std::vector<std::vector<Real>> temp_pos;
		//		temp_pos.emplace_back(his_ind[i]->variable().vect());
		//		if (transform_coordinate) {
		//			for (size_t j = 0; j < temp_pos.size(); ++j) {
		//				for (size_t k = 0; k < num_var; ++k) {
		//					temp_pos[j][k] -= t_origin[k];
		//				}
		//			}
		//		}
		//		if (shiftPlane(temp_pos, lower_weigh, 1)>0) {
		//			lower_inx.push_back(i);
		//		}
		//		if (shiftPlane(temp_pos, upper_weigh, 0)<0) {
		//			upper_inx.push_back(i);
		//		}
		//	}
		//	if (lower_inx.empty() || upper_inx.empty()) {
		//		islocal = false;
		//	}
		//	else {
		//		islocal = true;
		//	}
		//	//再从两边各自找出离前沿个体最近的个体
		//}
		//return islocal;

  //      // method 4: 先基于密度聚类，再判断各类与边界的距离
		//auto& present_ind = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		//size_t num_var = CAST_CONOP(pro)->numberVariables();
		//if (present_ind.size() < num_var) {
		//	islocal = false;
		//}
		//else {
		//	size_t minimum_num_cluster = 2;//每个类中最小的个体数
		//	Real epsilon = 0.5 * 0.5;
		//	std::vector<std::vector<Real>> cluster_sol;
		//	for (size_t i = 0; i < present_ind.size(); ++i) {
		//		cluster_sol.emplace_back(present_ind[i]->variable().vect());
		//	}
		//	DBSCAN dscluster(minimum_num_cluster, epsilon, cluster_sol);
		//	dscluster.run();
		//	std::vector<int> cluster_id;//每个解所属的类
		//	for (size_t i = 0; i < dscluster.m_points.size(); ++i) {
		//		cluster_id.push_back(dscluster.m_points[i]->clusterID);
		//	}
		//	std::vector<int> cluster_num;//所有的类
		//	cluster_num.push_back(cluster_id[0]);
		//	for (size_t i = 1; i < cluster_id.size(); ++i) {
		//		if (std::find(cluster_num.begin(), cluster_num.end(), cluster_id[i]) == cluster_num.end()) {
		//			cluster_num.push_back(cluster_id[i]);
		//		}
		//	}
		//	//将相同的类的点放在一起
		//	std::vector<std::vector<size_t>> ind2cluster;
		//	for (size_t i = 0; i < cluster_num.size(); ++i) {
		//		std::vector<size_t> temp_cluster;
		//		for (size_t j = 0; j < cluster_id.size(); ++j) {
		//			if (cluster_id[j] == cluster_num[i]) {
		//				temp_cluster.push_back(j);
		//			}
		//		}
		//		ind2cluster.emplace_back(temp_cluster);
		//	}

		//	//根据每个类前沿的位置判断
		//	//先得到边界
		//	std::vector<std::pair<Real, Real>> bound;
		//	for (size_t i = 0; i < num_var; ++i) {
		//		bound.emplace_back(CAST_CONOP(pro)->domain().range(i).limit);
		//	}
		//	bool flag = false;
		//	for (size_t i = 0; i < ind2cluster.size(); ++i) {
		//		auto inds_inx = ind2cluster[i];
		//		std::vector<Real> center_pos;
		//		for (size_t i = 0; i < num_var; ++i) {
		//			Real temp_sum = 0.;
		//			for (size_t j = 0; j < inds_inx.size(); ++j) {
		//				temp_sum += (present_ind[inds_inx[j]]->variable()[i]);
		//			}
		//			center_pos.push_back(temp_sum / inds_inx.size());
		//		}
		//		//中心点到bound的最短距离及维度
		//		flag = false;
		//		std::vector<Real> min_dim_dist(num_var, 0.);
		//		for (size_t j = 0; j < num_var; ++j) {
		//			min_dim_dist[j] = std::fabs(center_pos[j] - bound[j].first) < std::fabs(center_pos[j] - bound[j].second) ? std::fabs(center_pos[j] - bound[j].first) : std::fabs(center_pos[j] - bound[j].second);
		//			Real delta = min_dim_dist[j] / (bound[j].second - bound[j].first);
		//			if (delta < 0.1) {
		//				flag = true;
		//				break;
		//			}
		//		}
		//	}
		//	if (!flag) {
		//		islocal = true;
		//	}
		//}
		//return islocal;

	}

	//判断一个集合的点是否被一个超平面覆盖
	Real SPMOEA0::shiftPlane(std::vector<std::vector<Real>>& inputdata, std::vector<Real> weigh, bool flag) {
		Real offset = 0.;
		std::vector<Real> all_offset;
		for (size_t i = 0; i < inputdata.size(); ++i) {
			Real sum = 0.;
			for (size_t j = 0; j < weigh.size() - 1; ++j) {
				sum += (inputdata[i][j] * weigh[j]);
			}
			sum += weigh.back();
			sum -= inputdata[i].back();
			all_offset.push_back(sum);
		}
		if (flag) {//需要所有点在超平面下面，offset>0
			Real min_pos_off=(Real)INT16_MAX;
			Real min_neg_off= (Real)INT16_MAX;
			for (size_t i = 0; i < all_offset.size(); ++i) {
				if (all_offset[i] >= 0) {
					min_pos_off = min_pos_off < all_offset[i] ? min_pos_off : all_offset[i];
				}
				else {
					min_neg_off = min_neg_off < all_offset[i] ? min_neg_off : all_offset[i];
				}
			}
			if (min_neg_off<0) {
				flag = false;
				offset = min_neg_off-0.000001;
			}
			else {
				flag = true;
				offset = min_pos_off-0.000001;
			}
		}
		else {//需要所有点在超平面上面，offset<0
			Real max_pos_off = -1.*INT16_MAX;
			Real max_neg_off = -1*INT16_MAX;
			for (size_t i = 0; i < all_offset.size(); ++i) {
				if (all_offset[i] >= 0) {
					max_pos_off = max_pos_off >= all_offset[i] ? max_pos_off : all_offset[i];
				}
				else {
					max_neg_off = max_neg_off >= all_offset[i] ? max_neg_off : all_offset[i];
				}
			}
			if (max_pos_off >0) {
				flag = false;
				offset = max_pos_off+0.000001;
			}
			else {
				flag = true;
				offset = max_neg_off+0.000001;
			}
		}
		return offset;
	}

	void SPMOEA0::updateMultiPop(const std::vector<std::vector<size_t>>& clusters, const std::vector<size_t>& pre_front, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<size_t> front_space;
		for (size_t i = 0; i < clusters.size() - 1; ++i) {
			if (clusters[i].size() == 1) {
				front_space.push_back(clusters[i][0]);
			}
		}
		std::vector<size_t> delete_pop_inx;
		std::vector<size_t> preserve_pop_inx;
		for (size_t i = 0; i < pre_front.size(); ++i) {
			if (std::find(front_space.begin(), front_space.end(), pre_front[i]) == front_space.end()) {
				delete_pop_inx.push_back(pre_front[i]);
			}
			else {
				preserve_pop_inx.push_back(pre_front[i]);
			}
		}
		std::vector<size_t> add_pop_inx;
		for (size_t i = 0; i < front_space.size(); ++i) {
			if (std::find(pre_front.begin(), pre_front.end(), front_space[i]) == pre_front.end()) {
				add_pop_inx.push_back(front_space[i]);
			}
		}
		getPop().remove(getPop().end()-1);//删除最后的探索子种群
		//当前种群中没有该子种群，再添加
		for (size_t i = 0; i < add_pop_inx.size(); ++i) {
			std::vector<size_t> temp_space;
			temp_space.push_back(add_pop_inx[i]);

			bool pop_exist = false;
			for (size_t j = 0; j < getPop().size(); ++j) {
				auto temp = getPop()[j].getLocateSubspace();
				if (std::find(temp.begin(), temp.end(), add_pop_inx[i]) != temp.end()) {
					pop_exist = true;
					break;
				}
			}
			if (!pop_exist) {
				generateGoodPop(temp_space, pro, alg, rnd);
			}
		}
	}

	//bool SPMOEA0::setLinks(std::vector<std::vector<Real>>& v1, std::vector<std::vector<Real>>& v2) {
	//	//求两类点之间的最短距离
	//	Real min_dist=(Real)INT16_MAX;
	//	std::map<size_t,size_t> inx_links;
	//	size_t inx1;
	//	for (size_t i = 0; i < v1.size(); ++i) {
	//		Real temp_dist=(Real)INT16_MAX;
	//		size_t inx2;
	//		for (size_t j = 0; j < v2.size(); ++j) {
	//			auto dist = euclideanDistance(v1[i].begin(),v1[i].end(),v2[j].begin());
	//			if (dist < temp_dist) {
	//				temp_dist = dist;
	//				inx2 = j;
	//			}
	//		}
	//		inx_links.insert(std::make_pair(i,inx2));
	//		if (min_dist < temp_dist) {
	//			min_dist = temp_dist;
	//			inx1 = i;
	//		}
	//	}
	//	auto p1 = v1[inx1];
	//	auto p2 = v2[inx_links[inx1]];
	//	//每个集合的形心
	//}

	std::vector<std::pair<size_t, size_t>> SPMOEA0::frontChanged() {
		//先根据当前子空间信息找出前沿子空间
		std::vector<size_t> front_spaces_inx;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_spaces_inx.push_back(i);
			}
		}
		auto pre_front_clusters = getMO_HLC().getFrontClusters();
		//再与当前聚类子空间比较
		auto cur_front_clusters = getMO_HLC().clusterFrontSpace(front_spaces_inx,false);
		bool flag = false;
		std::vector<std::pair<size_t, size_t>> cluster_match;
		for (size_t i = 0; i < pre_front_clusters.size(); ++i) {
			for (size_t j = 0; j < cur_front_clusters.size(); ++j) {
				bool same = ifSame(pre_front_clusters[i],cur_front_clusters[j]);
				if (same) {
					cluster_match.emplace_back(std::make_pair<>(i,j));
				}
			}
		}
		return cluster_match;
	}

	std::vector<std::vector<size_t>> SPMOEA0::clusterSubspaces() {
		//找出前沿子空间和非前沿子空间
		std::vector<std::vector<size_t>> clusters;
		std::vector<size_t> front_space_inx;
		std::vector<size_t> behind_space_inx;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_space_inx.push_back(i);
			}
			else {
				behind_space_inx.push_back(i);
			}
		}
		for (size_t i = 0; i < front_space_inx.size(); ++i) {
			std::vector<size_t> temp;
			temp.push_back(front_space_inx[i]);
			clusters.emplace_back(temp);
		}
		//auto front_clusters = getMO_HLC().clusterFrontSpace(front_space_inx, false);
		//std::vector<std::vector<size_t>> clusters = front_clusters;
		
		getMO_HLC().setFrontClusters(clusters);
		clusters.emplace_back(behind_space_inx);
		//getMO_HLC().setSpaceClusters(clusters);
		return clusters;
	}

	std::vector<std::vector<size_t>> SPMOEA0::spaceCluster() {
		//找出前沿子空间和非前沿子空间
		std::vector<size_t> front_space_inx;
		std::vector<size_t> behind_space_inx;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_space_inx.push_back(i);
			}
			else {
				behind_space_inx.push_back(i);
			}
		}
		//生成多种群
		auto space_clusters = getMO_HLC().clusterFrontSpace(front_space_inx, true);
		auto front_clusters = getMO_HLC().clusterFrontSpace(front_space_inx, false);
		std::vector<size_t> last_cluster;
		for (size_t i = 0; i < behind_space_inx.size(); ++i) {
			bool flag = true;
			for (size_t j = 0; j < space_clusters.size(); ++j) {
				if (std::find(space_clusters[j].begin(), space_clusters[j].end(), behind_space_inx[i]) != space_clusters[j].end()) {
					flag = false;
					break;
				}
			}
			if (flag) {
				last_cluster.push_back(behind_space_inx[i]);
			}
		}
		space_clusters.emplace_back(last_cluster);
		getMO_HLC().setFrontClusters(front_clusters);
		//getMO_HLC().setSpaceClusters(space_clusters);
		return space_clusters;
	}

	std::vector<Real> SPMOEA0::localMutationSol(const std::vector<Real>& sol, const std::vector<std::pair<Real, Real>>& boundary, Problem *pro, Random *rnd) {
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
		//repairSol(inds, pro, rnd);
		return inds;
	}

	std::vector<Real> SPMOEA0::vectorMutationSol(const std::vector<Real>& sol1, const std::vector<Real>& sol2, Problem *pro, Random *rnd) {
		//先求穿出点坐标
		std::vector<Real> boundary_points;
		for (size_t i = 0; i < sol1.size(); ++i) {
			boundary_points.clear();
			Real temp_ratio1 = (m_pop_var_range[i].first - sol1[i]) / (sol1[i] - sol2[i]);
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
			temp_off.push_back(sol1[i] + temp_rand * (boundary_points[i] - sol1[i]));
		}
		//repairSol(temp_off, pro, rnd);
		return temp_off;
	}

	/*void SPMOEA0::repairSol(std::vector<Real>& sol, std::vector<size_t> space_inx,Random *rnd) {
		bool flag = false;
		for (size_t i = 0; i < space_inx.size(); ++i) {
			auto bound = getMO_HLC().subspaceTree().getBox(space_inx[i]);
			if (normalSol(sol, bound)) {
				flag=true;
				break;
			}
		}
		if (!flag) {
			size_t inx = (size_t)std::floor(rnd->uniform.next() * space_inx.size());
			auto bound= getMO_HLC().subspaceTree().getBox(space_inx[inx]);
			for (size_t i = 0; i < sol.size(); ++i) {
				sol[i] = bound[i].first + rnd->uniform.next() * (bound[i].second - bound[i].first);
			}
		}
	}*/

	void SPMOEA0::updateObjSpaceInfo(const Population<Solution<>>& pop, Problem *pro, bool changed) {
		/*更新目标值范围*/
		updateSubObjOpt(pop);
		bool b = updateObjRange(pop, pro);
		/*更新目标空间划分*/
		if (changed && !m_normalize) {
			updateObjTree();
		}
		/**********************************/
		/*      目标子空间加入个体信息    */
		/**********************************/
		///*目标空间清空,搜索空间选择性清空*/
		//m_obj_region_info.clear();
		//for (size_t i = 0; i < m_num_region_obj; ++i) {
		//	m_obj_region_info.emplace_back(new ObjRegionInfo);
		//}
		///*先更新目标子空间*/
		//for (size_t i = 0; i < pop.size(); i++) {
		//	auto temp_obj = pop[i].objective();
		//	if (m_normalize) {
		//		for (size_t j = 0; j < temp_obj.size(); ++j) {
		//			temp_obj[j] = (temp_obj[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
		//		}
		//	}
		//	size_t idx_obj_region = m_mo_hlc->getObjspaceTree().getRegionIdx(temp_obj);
		//	if (m_obj_region_info[idx_obj_region]->obj_optima.empty()) {//obj
		//		for (size_t j = 0; j < temp_obj.size(); j++) {
		//			m_obj_region_info[idx_obj_region]->obj_optima.emplace_back(pop[i].objective()[j]);
		//		}
		//	}
		//	else {
		//		for (size_t j = 0; j < temp_obj.size(); j++) {
		//			Real temp = m_obj_region_info[idx_obj_region]->obj_optima[j];
		//			m_obj_region_info[idx_obj_region]->obj_optima[j] = std::min(temp, pop[i].objective()[j]);
		//		}
		//	}
		//	m_obj_region_info[idx_obj_region]->m_curr_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));//solution
		//	//个体的索引
		//	m_obj_region_info[idx_obj_region]->ind_idx.push_back(i);
		//	if (pop[i].fitness() < m_obj_region_info[idx_obj_region]->m_obj_rank) { //目标子空间排序
		//		m_obj_region_info[idx_obj_region]->m_obj_rank = pop[i].fitness();
		//	}
		//}
	}

	void SPMOEA0::updateSubObjOpt(const Population<Solution<>>& pop) {
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

	bool SPMOEA0::updateObjRange(const Population<Solution<>>& pop, Problem *pro) {
		std::vector<std::pair<Real, Real>> pop_range_before = getPopRange();
		int M = CAST_CONOP(pro)->numberObjectives();
		int N = CAST_CONOP(pro)->numberVariables();
		//
		return false;
	}

	void SPMOEA0::updateObjTree() {
		auto obj_range = getFrontObjRange();
		getMO_HLC().getObjspaceTree().setInitBox(obj_range);
		size_t m_num_region_obj = numObjRegion();
		getMO_HLC().getObjspaceTree().inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		getMO_HLC().getObjspaceTree().buildIndex();
	}

	void SPMOEA0::updateSolSpaceInfo(const Population<Solution<>>& pop, Problem *pro, bool b) {
		//没有解和其他信息的子空间，能否由邻域信息推断？
		//得到优化模式
		//依次遍历解来更新子空间
		//有些信息需要及时清空
		if (b) {//b=true表示更新历史信息
			size_t num = pop.size() / 2;
			for (size_t i = 0; i < num; ++i) {
				auto temp_var = pop[i].variable().vect();
				auto temp_obj = pop[i].objective();
				size_t idx_var_region = getMO_HLC().subspaceTree().getRegionIdx(temp_var);
				//加入子目标最优值
				if (getMO_HLC().getSubspaceInfo(idx_var_region).m_history_inds.empty()) {//var
					for (size_t j = 0; j < temp_obj.size(); j++) {
						getMO_HLC().getSubspaceInfo(idx_var_region).m_subObj_optima.push_back(temp_obj[j]);
						//m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve.push_back(0.);
					}
				}
				else {
					for (size_t j = 0; j < temp_obj.size(); j++) {
						Real temp = getMO_HLC().getSubspaceInfo(idx_var_region).m_subObj_optima[j];
						if (temp > temp_obj[j]) {
							getMO_HLC().getSubspaceInfo(idx_var_region).m_subObj_optima[j] = temp_obj[j];
							//Real temp_improve = m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve[j];
							//m_mo_hlc->get_subspace_info()[idx_var_region]->m_subObj_improve[j] = std::min(temp_improve,temp_obj[j]-temp);
						}
					}
				}
				//update historical sols
				getMO_HLC().getSubspaceInfo(idx_var_region).m_history_inds.emplace_back(std::make_shared<Solution<>>(pop[i]));
				//update sample frequency
				getMO_HLC().getSubspaceInfo(idx_var_region).m_sub_freq++;
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
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				getMO_HLC().getSubspaceInfo(i).m_curr_sols.clear();
				//m_mo_hlc->getSubspaceInfo(i).m_curr_sols_idx.clear();
				getMO_HLC().getSubspaceInfo(i).m_best_rank = INT16_MAX;
			}
			int max_rank = 0;
			for (size_t i = 0; i < pop.size(); i++) {
				auto temp_var = pop[i].variable().vect();
				auto temp_obj = pop[i].objective();
				size_t idx_var_region = getMO_HLC().subspaceTree().getRegionIdx(temp_var);
				//加入子空间内的最好排序值
				if (getMO_HLC().getSubspaceInfo(idx_var_region).m_curr_sols.empty()) {//rank
					getMO_HLC().getSubspaceInfo(idx_var_region).m_best_rank = pop[i].fitness();
					if (pop[i].fitness() > max_rank) {
						max_rank = pop[i].fitness();
					}
				}
				else {
					getMO_HLC().getSubspaceInfo(idx_var_region).m_best_rank = std::min(getMO_HLC().getSubspaceInfo(idx_var_region).m_best_rank, pop[i].fitness());
					if (getMO_HLC().getSubspaceInfo(idx_var_region).m_best_rank > max_rank) {
						max_rank = getMO_HLC().getSubspaceInfo(idx_var_region).m_best_rank;
					}
				}
				//加入子空间内的当前解
				getMO_HLC().getSubspaceInfo(idx_var_region).m_curr_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));//solution
				//m_mo_hlc->getSubspaceInfo(idx_var_region).m_curr_sols_idx.push_back(i);
				//加入子空间内的前沿解
				if (pop[i].fitness() == 0) {
					getMO_HLC().getSubspaceInfo(idx_var_region).m_front_sol_in_subspace.push_back(std::make_shared<Solution<>>(pop[i]));
				}
			}

			//为没有当前个体的子空间分配rank值
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				auto& spaceinfo = getMO_HLC().getSubspaceInfo(i);
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
			calSpacePotential(obj_var_potential.first, pro);
			/*for (size_t i = 0; i < explore_space.size(); ++i) {
				m_mo_hlc->getSubspaceInfo(explore_space[i]).m_explore_potential = m_mo_hlc->calExploreSpacePotential(explore_space[i]);
			}
			for (size_t i = 0; i < exploit_space.size(); ++i) {
				m_mo_hlc->getSubspaceInfo(exploit_space[i]).m_exploit_potential = m_mo_hlc->calExploitSpacePotential(exploit_space[i], obj_var_potential);
			}*/
		}
	}

	void SPMOEA0::updateArchive(size_t num, Problem *pro) {
		//从历史前沿中选择一定数量的个体作为archive
		//三种选择方式：1、拥挤选择；2、参考向量选择；3、子空间跨度选择
		std::vector<size_t> select_index;
	}

	void SPMOEA0::updateHistorySol(Problem *pro) {
		/*std::vector<Solution<>> temp_pop;
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
		}*/
	}

	bool SPMOEA0::updateHisFrontObjRange() {
		/*if (m_history_front_sols.size() > 1) {
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
		}*/
		return false;
	}

	void SPMOEA0::calSpacePotential(std::vector<size_t>& spaces, Problem *pro) {
		//根据子空间的采样反馈系数、最佳排序值、子目标的最好排序值,以及目标空间的特征影响的潜力值
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			auto& spaceinfo = getMO_HLC().getSubspaceInfo(i);
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
			Real potential = 1. / std::pow(1 + nd_rank, 0.5) + 1. / std::pow(best_rank, 0.1);
			Real num = 30 * (CAST_CONOP(pro)->numberVariables() - 1);
			Real a = 10, b = 1, c = 10. / num, d = 7;
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

	void SPMOEA0::updateSubObjRank(Problem *pro) {
		//先得到这些子空间目标值的最差值
		std::vector<Real> worst_value;
		//然后为各子空间排序
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		std::vector<size_t> has_counted(obj_num);
		for (size_t i = 0; i < obj_num; ++i) {
			size_t count_sort = 0;
			for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
				Real temp_obj;
				getMO_HLC().getSubspaceInfo(j).m_subObj_optima_rank.resize(obj_num);
				if (getMO_HLC().getSubspaceInfo(j).m_subObj_optima.empty()) {
					continue;
				}
				else {
					temp_obj = getMO_HLC().getSubspaceInfo(j).m_subObj_optima[i];
					count_sort++;
				}
				size_t count = 1;
				for (size_t k = 0; k < getMO_HLC().numSubspace(); ++k) {
					if (k != j) {
						if (getMO_HLC().getSubspaceInfo(k).m_subObj_optima.empty()) {
							continue;
						}
						else if (temp_obj > getMO_HLC().getSubspaceInfo(k).m_subObj_optima[i]) {
							count++;
						}
					}
				}
				getMO_HLC().getSubspaceInfo(j).m_subObj_optima_rank[i] = count;
			}
			has_counted[i] = count_sort;
		}
		for (size_t i = 0; i < obj_num; ++i) {
			for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
				if (getMO_HLC().getSubspaceInfo(j).m_subObj_optima.empty()) {
					getMO_HLC().getSubspaceInfo(j).m_subObj_optima_rank[i] = has_counted[i] + 1;
				}
			}
		}
	}

	bool SPMOEA0::ifApproxiConverge(Random *rnd) {
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

	std::vector<std::vector<size_t>> SPMOEA0::clusterExploitSpace(const KDTree& tree, const std::map<size_t, std::vector<size_t>>& frontspace) {
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
		std::vector<size_t> select_flag(front_space_inx.size(), 0);
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

	std::vector<size_t> SPMOEA0::getExploitSpace() {
		std::vector<size_t> exploit_space;
		//先找出子空间排序值为0的子空间
		std::vector<size_t> front_space;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_space.push_back(i);
			}
		}
		for (size_t i = 0; i < front_space.size(); ++i) {
			exploit_space.push_back(front_space[i]);
		}
		if (m_add_neighbor) {
			//再看其他子空间是否为这些子空间的邻域
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				for (size_t j = 0; j < front_space.size(); ++j) {
					std::list<size_t> neighbors;
					getMO_HLC().subspaceTree().findNeighbor(front_space[j], neighbors);
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

	std::vector<size_t> SPMOEA0::getExploreSpace() {
		//先找出子空间排序值为0的子空间
		std::vector<size_t> front_space;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_space.push_back(i);
			}
		}
		//先将非前沿子空间加入
		std::vector<size_t> explore_space;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (std::find(front_space.begin(), front_space.end(), i) == front_space.end()) {
				explore_space.push_back(i);
			}
		}

		if (m_add_neighbor) {
			//将已加入的子空间中的前排空间的邻域去掉
			std::vector<size_t> select_flag(explore_space.size(), 0);
			for (size_t i = 0; i < explore_space.size(); ++i) {
				for (size_t j = 0; j < front_space.size(); ++j) {
					std::list<size_t> neighbors;
					getMO_HLC().subspaceTree().findNeighbor(front_space[j], neighbors);
					if (std::find(neighbors.begin(), neighbors.end(), explore_space[i]) == neighbors.end()) {
						if (j == front_space.size() - 1) {
							select_flag[i] = 1;
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

	std::vector<std::vector<size_t>> SPMOEA0::clusterExploreSpace(const KDTree& tree) {
		//先得到探索子空间索引
		std::vector<size_t> explore_spaces = getExploreSpace();
		std::vector<std::vector<size_t>> explore_cluster;
		//再根据探索潜力值聚类
		//explore_cluster = getMO_HLC().clusterSubspace(explore_spaces,tree.size(), m_add_neighbor);
		//直接作为一类
		explore_cluster.push_back(explore_spaces);
		return explore_cluster;
	}

	std::vector<size_t> SPMOEA0::assignBasinResource(const std::vector<std::vector<size_t>>& potential_inds) {
		////归一化的潜力值
		std::vector<Real> resource;
		std::vector<size_t> assign_resource;
		//size_t total_num = 0;
		//Real sum = 0.;
		///*for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
		//	auto basin = m_mo_hlc->getBasinInfo(i);
		//	if (basin.flag == "explore") {
		//		resource.push_back(basin.m_basin_explore_potential);
		//	}
		//	else {
		//		resource.push_back(basin.m_basin_exploit_potential);
		//	}
		//	sum += resource.back();
		//}*/
		//for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
		//	auto basin = m_mo_hlc->getBasinInfo(i);
		//	resource.push_back(basin.m_basin_potential);
		//	sum += resource.back();
		//}
		//for (size_t i = 0; i < potential_inds.size(); ++i) {
		//	Real temp_potential = 0.1;
		//	resource.push_back(temp_potential);
		//	sum += temp_potential;
		//}
		//for (auto& i : resource) {
		//	i = i / sum;
		//}
		//for (size_t i = 0; i < resource.size(); ++i) {
		//	auto tmp = std::ceil(m_pop->size() * resource[i]);
		//	assign_resource.push_back(tmp);
		//	total_num += tmp;
		//}
		////去除多余个体
		//while (total_num > m_pop->size()) {
		//	size_t max = 0;
		//	size_t inx = 0;
		//	for (size_t i = 0; i < resource.size(); ++i) {
		//		if (assign_resource[i] > max) {
		//			max = assign_resource[i];
		//			inx = i;
		//		}
		//	}
		//	assign_resource[inx] = assign_resource[inx] - 1;
		//	total_num -= 1;
		//}
		return assign_resource;
	}

	bool SPMOEA0::popConverged() {
		//size_t count = 0;
		//for (size_t i = 0; i < m_pop->getOffspring().size(); i++) {
		//	if (m_pop->getOffspring()[i].fitness() == 0)//连续存活n代
		//		count++;
		//}
		//if (count > m_pop->size())
		//	return true;
		//else
		//	return false;
		return false;
	}

	void SPMOEA0::initiObjSpace(Problem *pro) {
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
		//m_mo_hlc->getObjspaceTree().setInitBox(m_obj_boundary);
		//m_mo_hlc->getObjspaceTree().inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		//m_mo_hlc->getObjspaceTree().buildIndex();
		///*m_obj_tree->setInitBox(m_obj_boundary);
		//m_obj_tree->inputRatioData(std::vector<Real>(m_num_region_obj, 1.0 / m_num_region_obj));
		//m_obj_tree->buildIndex();*/

		//updateObjSpaceInfo(m_pop->getParent(), pro, false);
		//std::vector<size_t> space_idx;
		//for (size_t i = 0; i < numObjSubspace(); ++i) {
		//	auto& spaceinfo = getObjRegionInfo(i);
		//	if (spaceinfo.m_obj_rank == 0) {
		//		space_idx.push_back(i);
		//	}
		//}
		//m_front_space_inx.push_back(space_idx);
	}

	//void SPMOEA0::initiVarSpace(Problem *pro) {
	//	/*************************************/
	//	/*      MO_HLC搜索空间划分初始化     */
	//	/*************************************/
	//	MOSP::initiVarSpace(pro);
	//}

	//std::vector<size_t> SPMOEA0::predictObjSpacePotential(Problem *pro) {//输出的是决策子空间的位置
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
	//			Real dist = (Real)INT_MAX;
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

	std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> SPMOEA0::predictObjSpacePotential(Problem *pro) {//输出的是决策子空间的位置
		//在目标空间上表现为：与已收敛子空间处于互不支配地位的子空间且里面个体很少，还有局部PF的边界点周围
		//在搜索空间上表现为：目标空间上述点在搜索空间对应点的邻域范围内
		std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> output_pair;
		//std::vector<std::vector<size_t>> potential_ind_idx;//第一行为所有目标的极值，剩余的为各子目标下的稀疏点对
		//std::vector<size_t> first_indi_index;//第一排个体索引
		//for (size_t i = 0; i < m_pop->size(); i++) {
		//	if (m_pop->at(i).fitness() == 0)
		//		first_indi_index.emplace_back(i);
		//}
		//size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		//std::vector<size_t> obj_point_index;
		////先加入前排个体中子目标最差的边界值索引
		//for (size_t i = 0; i < num_obj; i++) {
		//	std::vector<size_t> boundary_idx;
		//	size_t temp = first_indi_index[0];
		//	for (size_t j = 0; j < first_indi_index.size(); j++) {
		//		if (m_pop->at(temp).objective()[i] < m_pop->at(first_indi_index[j]).objective()[i]) {
		//			temp = first_indi_index[j];
		//		}
		//	}
		//	boundary_idx.push_back(temp);

		//	if (obj_point_index.empty() || std::find(obj_point_index.begin(), obj_point_index.end(), temp) == obj_point_index.end()) {
		//		obj_point_index.emplace_back(temp);
		//	}
		//	potential_ind_idx.push_back(boundary_idx);
		//}

		////再加入第一排稀疏个体的索引，计算距离前需要对目标空间进行归一化
		//// 在每个子目标上检索距离最远的个体
		//std::vector<size_t> sparse_idx;
		//std::vector<std::vector<Real>> temp_first_obj;
		//for (size_t i = 0; i < first_indi_index.size(); ++i) {
		//	auto temp_obj = m_pop->at(first_indi_index[i]).objective();
		//	if (m_normalize) {
		//		for (size_t j = 0; j < temp_obj.size(); ++j) {
		//			temp_obj[j] = (temp_obj[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
		//		}
		//	}
		//	temp_first_obj.push_back(temp_obj);
		//}
		//if (first_indi_index.size() > 1) {
		//	for (size_t i = 0; i < num_obj; ++i) {
		//		//所有解按照某一维目标排序
		//		//先构造元组
		//		std::vector<size_t> sparse_pair_idx;//存放每个目标上的稀疏点对
		//		std::vector<std::tuple<Real, size_t>> temp_dim_obj;
		//		for (size_t j = 0; j < temp_first_obj.size(); ++j) {
		//			temp_dim_obj.emplace_back(std::make_tuple<>(temp_first_obj[j][i], first_indi_index[j]));
		//		}
		//		//按目标值递增排序
		//		for (size_t j = 0; j < temp_dim_obj.size(); ++j) {
		//			for (size_t k = j + 1; k < temp_dim_obj.size(); ++k) {
		//				if (std::get<0>(temp_dim_obj[j]) > std::get<0>(temp_dim_obj[k])) {
		//					auto temp = temp_dim_obj[j];
		//					temp_dim_obj[j] = temp_dim_obj[k];
		//					temp_dim_obj[k] = temp;
		//				}
		//			}
		//		}
		//		//求最大间隔
		//		Real max_dist = 0;
		//		size_t idx = 0;
		//		for (size_t j = 1; j < temp_dim_obj.size(); ++j) {
		//			auto temp = std::get<0>(temp_dim_obj[j]) - std::get<0>(temp_dim_obj[j - 1]);
		//			if (temp > max_dist) {
		//				max_dist = temp;
		//				idx = j;
		//			}
		//		}
		//		sparse_pair_idx.push_back(std::get<1>(temp_dim_obj[idx - 1]));
		//		sparse_pair_idx.push_back(std::get<1>(temp_dim_obj[idx]));
		//		//判断是否已有相同的点对加入
		//		if (potential_ind_idx.size() == num_obj) {
		//			potential_ind_idx.push_back(sparse_pair_idx);
		//		}
		//		else {
		//			bool flag = false;
		//			for (size_t j = num_obj; j < potential_ind_idx.size(); ++j) {
		//				auto obj = potential_ind_idx[j];
		//				if (obj[0] == sparse_pair_idx[0] && obj[1] == sparse_pair_idx[1] || \
		//					obj[0] == sparse_pair_idx[1] && obj[1] == sparse_pair_idx[0]) {
		//					flag = true;
		//					break;
		//				}
		//			}
		//			if (!flag) {
		//				potential_ind_idx.push_back(sparse_pair_idx);
		//			}
		//		}

		//		if (sparse_idx.empty() || std::find(sparse_idx.begin(), sparse_idx.end(), std::get<1>(temp_dim_obj[idx])) == sparse_idx.end()) {
		//			sparse_idx.push_back(std::get<1>(temp_dim_obj[idx]));
		//		}
		//		if (sparse_idx.empty() || std::find(sparse_idx.begin(), sparse_idx.end(), std::get<1>(temp_dim_obj[idx - 1])) == sparse_idx.end()) {
		//			sparse_idx.push_back(std::get<1>(temp_dim_obj[idx - 1]));
		//		}
		//	}
		//}
		//for (size_t i = 0; i < sparse_idx.size(); ++i) {
		//	if (std::find(obj_point_index.begin(), obj_point_index.end(), sparse_idx[i]) == obj_point_index.end()) {
		//		obj_point_index.push_back(sparse_idx[i]);
		//	}
		//}

		////最后根据个体索引得到搜索子空间索引
		//std::vector<size_t> var_space_idx;
		//for (size_t i = 0; i < obj_point_index.size(); ++i) {
		//	auto position = m_pop->at(obj_point_index[i]).variable().vect();
		//	auto index = m_mo_hlc->subspaceTree().getRegionIdx(position);
		//	if (var_space_idx.empty() || std::find(var_space_idx.begin(), var_space_idx.end(), index) == var_space_idx.end()) {
		//		var_space_idx.push_back(index);
		//	}
		//}
		//output_pair.first = var_space_idx;
		//output_pair.second = potential_ind_idx;
		return output_pair;
	}

	//void SPMOEA0::updateVarSpacePotential() {
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

	/*std::vector<size_t> SPMOEA0::predictVarSpacePotential(Problem *pro) {
		std::vector<size_t> space_potential;

		return space_potential;
	}*/

	void SPMOEA0::envirSelection(Problem *pro, Random *rnd) {//是否趋近收敛有不同的环境选择策略
		////m_pop->envirSelection(pro, rnd, b);
		////基于子空间跨度的多参考选择方法
		//if (m_subspace_select) {
		//	auto select_pop = selectIndi(m_pop->getOffspring(), m_pop->size(), pro, rnd);
		//	for (size_t i = 0; i < m_pop->size(); ++i) {
		//		m_pop->at(i) = m_pop->getOffspring()[select_pop[i]];
		//		m_pop->at(i).setCounter(m_pop->at(i).surviveAge() + 1);
		//	}
		//}
		//else {
		//	m_pop->envirSelection(pro, rnd, false);
		//}
	}

	//void SPMOEA0_pop::envirSelection(Problem *pro, Random *rnd, bool b) {
	//	if (b) {//趋近收敛时的环境选择基于子空间
	//		//基于目标子空间进行选择，选择依据为：前沿子空间的均匀性和支配空间的收敛性
	//		//先对父代和子代进行非支配排序
	//		/*for (size_t i = 0; i < this->m_individuals.size(); ++i) {
	//			m_offspring[i + this->m_individuals.size()] = *this->m_individuals[i];
	//		}*/
	//		nondominatedSorting(m_offspring);
	//		/*auto select_pop = selectIndi(m_offspring, this->m_individuals.size(),pro);
	//		for (size_t i = 0; i < this->m_individuals.size(); ++i) {
	//			*this->m_individuals[i] = m_offspring[select_pop[i]];
	//			this->m_individuals[i]->setCounter(this->m_individuals[i]->surviveAge() + 1);
	//		}*/
	//	}
	//	else {
	//		survivorSelection(*this, m_offspring);
	//		for (size_t i = 0; i < this->m_individuals.size(); ++i) {
	//			this->m_individuals[i]->setCounter(this->m_individuals[i]->surviveAge() + 1);
	//		}
	//	}
	//}

	//std::vector<Real> SPMOEA0_pop::optimaPredict(const std::vector<std::vector<Real>>& points) {//预测得到的是坐标,points为一个点序列
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

	//std::vector<std::vector<std::vector<Real>>> SPMOEA0_pop::indi_time_series(const std::vector<std::vector<Solution<>>>& p) {//对memeory中个体直接找距离最小值形成序列
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

	//std::vector<size_t> SPMOEA0::selectIndi(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd) {
	//	size_t num_obj = CAST_CONOP(pro)->numberObjectives();
	//	auto pop_att = popAttach(pop);
	//	std::vector<size_t> select_index;//已经选择个体的索引
	//	std::vector<size_t> behind_space_index;//后排子空间索引
	//	std::vector<size_t> front_space_index;//前排子空间索引
	//	size_t front_space_num = 0;//前沿子空间内的个体总数
	//	std::vector<size_t> first_ind_index;//先得到第一排个体的索引
	//	for (size_t i = 0; i < pop.size(); ++i) {
	//		if (pop[i].fitness() == 0)
	//			first_ind_index.emplace_back(i);
	//	}
	//	for (const auto& ind : pop_att[0][0]) {
	//		front_space_index.emplace_back(ind.first);
	//		front_space_num += ind.second.size();
	//	}
	//	for (const auto& ind : pop_att[0][1]) {
	//		behind_space_index.emplace_back(ind.first);
	//	}
	//	std::map<size_t, std::vector<size_t>> map_first_indi;//得到前排个体分布在哪些子空间
	//	for (const auto& sp : pop_att[0][0]) {
	//		std::pair<size_t, std::vector<size_t>> temp;
	//		temp.first = sp.first;
	//		for (const auto& in : sp.second) {
	//			if (pop[in].fitness() == 0) {
	//				temp.second.push_back(in);
	//			}
	//		}
	//		map_first_indi.insert(temp);
	//	}
	//	std::map<size_t, std::vector<size_t>> map_selected_indi;//前沿子空间已经选的解，顺序索引
	//	std::map<size_t, std::vector<size_t>> map_no_selected_indi;//前沿子空间已经选的解，顺序索引
	//	//根据前排子空间内部的个体总数进行个体选择，选择时考虑个体年龄
	//	if (first_ind_index.size() > select_num) {//前沿个体充足时，进行子空间多参考点选择
	//		//先得到当前子目标的极值解的索引
	//		std::vector<size_t> sub_opt_idx;
	//		for (size_t j = 0; j < num_obj; ++j) {
	//			std::vector<Real> obj_value;
	//			for (size_t k = 0; k < first_ind_index.size(); ++k) {
	//				obj_value.push_back(pop[first_ind_index[k]].objective()[j]);
	//			}
	//			size_t min_idx = min_element(obj_value.begin(), obj_value.end()) - obj_value.begin();
	//			sub_opt_idx.push_back(min_idx);
	//		}
	//		//再选择子空间内某个子目标的极值点
	//		for (const auto& i : map_first_indi) {
	//			std::pair<size_t, std::vector<size_t>> temp1;
	//			std::pair<size_t, std::vector<size_t>> temp2;
	//			temp1.first = i.first;
	//			temp2.first = i.first;
	//			auto ind_idx = i.second;
	//			if (ind_idx.size() == 1) {
	//				select_index.push_back(ind_idx[0]);
	//				temp1.second.push_back(ind_idx[0]);
	//			}
	//			else {//比较子空间内点子目标值
	//				std::vector<Real> obj_value;
	//				for (size_t k = 0; k < ind_idx.size(); ++k) {
	//					obj_value.push_back(pop[ind_idx[k]].objective()[0]);
	//				}
	//				size_t min_idx = min_element(obj_value.begin(), obj_value.end()) - obj_value.begin();
	//				select_index.push_back(ind_idx[min_idx]);
	//				temp1.second.push_back(min_idx);
	//			}
	//			std::vector<size_t> no_select_idx;
	//			for (size_t j = 0; j < map_first_indi[i.first].size(); ++j) {
	//				if (std::find(temp1.second.begin(), temp1.second.end(), j) == temp1.second.end()) {
	//					no_select_idx.push_back(j);
	//				}
	//			}
	//			temp2.second = no_select_idx;
	//			map_selected_indi.insert(temp1);
	//			map_selected_indi.insert(temp2);
	//		}
	//		std::map<size_t, std::vector<std::vector<Real>>> map_subspace_dist_matrix;//子空间未选点与其邻域已选点的距离矩阵
	//		std::map<size_t, std::vector<size_t>> map_neigh_select_ind;//邻域已选个体索引,真实索引
	//		std::map<size_t, std::vector<size_t>> map_neigh_inx;//邻域索引
	//		for (const auto& i : map_first_indi) {
	//			std::pair<size_t, std::vector<std::vector<Real>>> temp;
	//			std::pair<size_t, std::vector<size_t>> neigh_select_ind;
	//			std::pair<size_t, std::vector<size_t>> neigh_inx;
	//			neigh_select_ind.first = i.first;
	//			temp.first = i.first;
	//			neigh_inx.first = i.first;
	//			auto idx = i.second;
	//			//先找出邻域子空间索引，并得到所有邻域已选个体索引
	//			std::list<size_t> neighbors;
	//			std::vector<size_t> front_nei_select_idx;//邻域已选个体，实际索引
	//			std::vector<size_t> front_nei_idx;
	//			getMO_HLC().getObjspaceTree().findNeighbor(i.first, neighbors);
	//			for (const auto& ss : map_first_indi) {
	//				if (find(neighbors.begin(), neighbors.end(), ss.first) != neighbors.end()) {
	//					front_nei_idx.push_back(ss.first);
	//					//front_nei_select_idx.push_back(ss.first);
	//					auto inx = ss.second;
	//					auto ind_selected = map_selected_indi[ss.first];
	//					for (size_t j = 0; j < ind_selected.size(); ++j) {
	//						front_nei_select_idx.push_back(ss.second[ind_selected[j]]);
	//					}
	//				}
	//			}
	//			neigh_inx.second = front_nei_idx;
	//			map_neigh_inx.insert(neigh_inx);
	//			//子空间未选个体索引
	//			std::vector<size_t> no_select_inx;
	//			auto no_se = map_no_selected_indi[i.first];
	//			for (size_t j = 0; j < no_se.size(); ++j) {
	//				no_select_inx.push_back(i.second[no_se[j]]);
	//			}
	//			//未选点与已选点的距离矩阵，如果距离矩阵为空怎么办？
	//			std::vector<std::vector<Real>> dist_matrix;
	//			for (size_t j = 0; j < no_select_inx.size(); ++j) {
	//				std::vector<Real> ind_dist;
	//				for (size_t k = 0; k < front_nei_select_idx.size(); ++k) {
	//					ind_dist.push_back(euclideanDistance(pop[no_select_inx[j]].objective().begin(), pop[no_select_inx[j]].objective().end(), pop[front_nei_select_idx[k]].objective().begin()));
	//				}
	//				dist_matrix.push_back(ind_dist);
	//			}
	//			temp.second = dist_matrix;
	//			neigh_select_ind.second = front_nei_select_idx;
	//			map_subspace_dist_matrix.insert(temp);
	//			map_neigh_select_ind.insert(neigh_select_ind);
	//		}
	//		//然后根据每个子空间内个体最小距离的统计指标选择要选的子空间
	//		size_t select_subspace_idx = 0;
	//		while (select_index.size() < select_num) {
	//			//选择子空间内未选个体最小距离最大的子空间，从里面采用最小距离最大规则选出一个解
	//			std::vector<std::pair<size_t, Real>> sub_min_max;
	//			std::map<size_t, size_t> sub_min_max_inx;//每个子空间最大值的索引
	//			for (const auto& i : map_subspace_dist_matrix) {
	//				auto matri = i.second;
	//				Real max_dist = 0.;
	//				size_t inx = 0;
	//				if (matri.size() != 0) {
	//					std::vector<size_t> temp_min;
	//					for (size_t j = 0; j < i.second.size(); ++j) {
	//						auto min_temp_dist = *min_element(i.second[j].begin(), i.second[j].end());
	//						temp_min.push_back(min_temp_dist);
	//					}
	//					max_dist = *max_element(temp_min.begin(), temp_min.end());
	//					inx = max_element(temp_min.begin(), temp_min.end()) - temp_min.begin();
	//				}
	//				else {
	//					size_t a = 1;
	//				}
	//				std::pair<size_t, Real> temp1;
	//				std::pair<size_t, size_t> temp2;
	//				temp1.first = i.first;
	//				temp1.second = max_dist;
	//				temp2.first = i.first;
	//				temp2.second = inx;
	//				sub_min_max.push_back(temp1);
	//				sub_min_max_inx.insert(temp2);
	//			}
	//			//选出子空间
	//			size_t subspace_idx;
	//			std::vector<std::pair<size_t, Real>> subspace_dist = sub_min_max;
	//			for (size_t j = 0; j < sub_min_max.size() - 1; ++j) {
	//				for (size_t k = j + 1; k < sub_min_max.size(); ++k) {
	//					if (subspace_dist[j].second < subspace_dist[k].second) {
	//						std::pair<size_t, Real> temp_pair;
	//						temp_pair.first = subspace_dist[k].first;
	//						temp_pair.second = subspace_dist[k].second;
	//						subspace_dist[k].first = subspace_dist[j].first;
	//						subspace_dist[k].second = subspace_dist[j].second;
	//						subspace_dist[j].first = temp_pair.first;
	//						subspace_dist[j].second = temp_pair.second;
	//					}
	//				}
	//			}
	//			for (size_t i = 0; i < subspace_dist.size(); ++i) {
	//				if (map_selected_indi[subspace_dist[i].first].size() < map_first_indi[subspace_dist[i].first].size()) {
	//					subspace_idx = subspace_dist[i].first;
	//					break;
	//				}
	//			}
	//			size_t max_idx = sub_min_max_inx[subspace_idx];//指示未选点的索引
	//			auto temp_idx = map_no_selected_indi[subspace_idx][max_idx];
	//			//所选点的索引
	//			size_t sele_idx = map_first_indi[subspace_idx][temp_idx];//未选点的实际索引
	//			//更新邻域子空间所选点和距离矩阵
	//			auto nei_idx = map_neigh_inx[subspace_idx];
	//			for (const auto& i : map_first_indi) {
	//				if (std::find(nei_idx.begin(), nei_idx.end(), i.first) != nei_idx.end()) {
	//					if (i.first == subspace_idx) {
	//						//去除选中点所在行
	//						map_subspace_dist_matrix[i.first].erase(map_subspace_dist_matrix[i.first].begin() + max_idx);
	//						map_no_selected_indi[i.first].erase(map_no_selected_indi[i.first].begin() + max_idx);
	//						map_selected_indi[i.first].push_back(temp_idx);
	//					}
	//					//添加已选点的列
	//					std::vector<size_t> no_sele_ind;
	//					for (size_t j = 0; j < map_no_selected_indi[i.first].size(); ++j) {
	//						no_sele_ind.push_back(i.second[map_no_selected_indi[i.first][j]]);
	//					}
	//					for (size_t j = 0; j < map_subspace_dist_matrix[i.first].size(); ++j) {
	//						map_subspace_dist_matrix[i.first][j].push_back(euclideanDistance(pop[no_sele_ind[j]].objective().begin(), pop[no_sele_ind[j]].objective().end(), pop[sele_idx].objective().begin()));
	//					}
	//					map_neigh_select_ind[i.first].push_back(sele_idx);
	//				}
	//			}
	//			select_index.push_back(sele_idx);
	//		}
	//	}
	//	else if (front_space_num > select_num) {//前沿个体不足时，前沿子空间个体充足时
	//		bool flag = 1;
	//		if (flag) {//目标空间选择
	//			//先选前沿个体
	//			select_index = first_ind_index;
	//			//再依次选择前排子空间中的其他个体
	//			auto front_subspace_ind = pop_att[0][0];
	//			for (const auto& i : front_subspace_ind) {
	//				std::pair<size_t, std::vector<size_t>> temp;
	//				temp.first = i.first;
	//				std::vector<size_t> f_ind_idx;
	//				for (size_t j = 0; j < i.second.size(); ++j) {
	//					if (find(first_ind_index.begin(), first_ind_index.end(), i.second[j]) != first_ind_index.end()) {
	//						f_ind_idx.push_back(i.second[j]);
	//					}
	//				}
	//				temp.second = f_ind_idx;
	//				map_selected_indi.insert(temp);//前沿子空间已经选的解
	//			}
	//			while (select_index.size() < select_num) {
	//				for (const auto& i : front_subspace_ind) {
	//					if (map_selected_indi[i.first].size() < i.second.size()) {
	//						//从剩余个体中选择一个
	//						for (size_t j = 0; j < i.second.size(); ++j) {
	//							if (find(map_selected_indi[i.first].begin(), map_selected_indi[i.first].end(), i.second[j]) == map_selected_indi[i.first].end()) {
	//								select_index.push_back(i.second[j]);
	//								map_selected_indi[i.first].push_back(i.second[j]);
	//								break;
	//							}
	//						}
	//					}
	//					if (select_index.size() >= select_num) {
	//						break;
	//					}
	//				}
	//			}
	//		}
	//		else {//解空间选择，根据解空间分布选择，或者根据子空间的潜力值选择
	//			//将前沿子空间的点分成两部分，前沿解和非前沿解

	//			//将所有解映射到搜索空间，找出离前沿解最远的个体
	//		}
	//	}
	//	else {//前排个体数不足时，前排全选，后面的轮流选
	//		//先加入前排子空间个体的索引
	//		for (const auto& ind : pop_att[0][0]) {
	//			for (const auto& p : ind.second) {
	//				select_index.push_back(p);
	//			}
	//		}
	//		//再在后排子空间轮流选择
	//		std::vector<size_t> space_select_num(pop_att[0][1].size(), 0);//每个后排子空间已选择的个数
	//		std::vector<std::vector<size_t>> behind_spaces;
	//		for (const auto& sp : pop_att[0][1]) {
	//			behind_spaces.push_back(sp.second);
	//		}
	//		while (select_index.size() < select_num) {
	//			for (size_t i = 0; i < behind_spaces.size(); ++i) {
	//				if (space_select_num[i] < behind_spaces.size()) {
	//					size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
	//					if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
	//						select_index.push_back(behind_spaces[i][select_idx]);
	//						space_select_num[i]++;
	//					}
	//				}
	//				if (select_index.size() >= select_num) {
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	//if (front_space_num > m_pop->size()) {//前排个体充足时，进行子空间多参考点选择
	//	//	//计算前排子空间跨度
	//	//	std::map<size_t, size_t> space_select_num;
	//	//	std::vector<Real> subspace_span;
	//	//	for (const auto& subspace : pop_att[0][0]) {
	//	//		std::vector<Real> ideal_point;
	//	//		std::vector<Real> nadir_point;
	//	//		auto sub_boundary = m_obj_tree->getBox(subspace.first);
	//	//		for (size_t j = 0; j < sub_boundary.size(); ++j) {
	//	//			ideal_point.push_back(sub_boundary[j].first);
	//	//			nadir_point.push_back(sub_boundary[j].second);
	//	//		}
	//	//		std::vector<size_t> ind_index = subspace.second;
	//	//		std::vector<std::vector<Real>> sub_points;
	//	//		for (size_t q = 0; q < ind_index.size(); ++q) {
	//	//			auto temp = pop[ind_index[q]].objective();
	//	//			if (m_normalize) {
	//	//				for (size_t j = 0; j < temp.size(); ++j) {
	//	//					temp[j] = (temp[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
	//	//				}
	//	//			}
	//	//			sub_points.push_back(temp);
	//	//		}
	//	//		subspace_span.push_back(calSubspaceSpan(ideal_point, nadir_point, sub_points));
	//	//	}
	//	//	//为前排子空间分配个体数
	//	//	std::vector<size_t> subspace_ind_num;
	//	//	Real sum_span = 0.;
	//	//	for (size_t p = 0; p < subspace_span.size(); ++p) {
	//	//		sum_span += subspace_span[p];
	//	//	}
	//	//	for (size_t p = 0; p < subspace_span.size(); ++p) {
	//	//		subspace_ind_num.push_back(1 + std::ceil(m_pop->size() * subspace_span[p] / sum_span));
	//	//	}
	//	//	//从前排子空间选择个体，满足优先级，即先子目标极值，再前排个体，最后其他个体
	//	//	size_t index = 0;
	//	//	m_over_space_index.clear();
	//	//	for (const auto& subspace : map_first_indi) {
	//	//		std::vector<Real> anchor_point;
	//	//		auto sub_boundary = m_obj_tree->getBox(subspace.first);
	//	//		for (size_t p = 0; p < sub_boundary.size(); ++p) {
	//	//			anchor_point.push_back(sub_boundary[p].first);
	//	//		}
	//	//		std::vector<size_t> ind_index = subspace.second;
	//	//		/*std::vector<std::vector<Real>> sub_points;
	//	//		for (size_t p = 0; p < ind_index.size(); ++p) {
	//	//			sub_points.push_back(pop[ind_index[p]].objective());
	//	//		}*/
	//	//		auto select_ind=selectIndSubspace(anchor_point, ind_index, subspace_ind_num[index]);
	//	//		std::vector<size_t> select_ind_index = std::get<0>(select_ind);
	//	//		bool flag= std::get<1>(select_ind);
	//	//		if (flag) {
	//	//			m_over_space_index.push_back(subspace.first);
	//	//		}
	//	//		for (size_t p = 0; p < select_ind_index.size(); ++p) {
	//	//			select_index.push_back(select_ind_index[p]);
	//	//		}
	//	//		index++;
	//	//		space_select_num.insert(std::make_pair<>(subspace.first,select_ind_index.size()));
	//	//	}
	//	//	//当选择个体超出时，进行拥挤淘汰，1216行越界
	//	//	while (select_index.size() > m_pop->size())
	//	//	{
	//	//		//先找出子目标极值个体
	//	//		std::vector<size_t> sub_optima_inx;
	//	//		for (size_t p = 0; p < num_obj; ++p) {
	//	//			size_t inx = 0;
	//	//			Real min_v = 1. * 10e14;
	//	//			for (size_t q = 0; q < select_index.size(); ++q) {
	//	//				if (pop[select_index[q]].objective()[p] < min_v) {
	//	//					min_v = pop[select_index[q]].objective()[p];
	//	//					inx = q;
	//	//				}
	//	//			}
	//	//			if (sub_optima_inx.empty() || std::find(select_index.begin(), select_index.end(), inx) == select_index.end()) {
	//	//				sub_optima_inx.push_back(select_index[inx]);
	//	//			}
	//	//		}
	//	//		//再计算个体的拥挤度，只计算与邻域个体的距离，有问题
	//	//		std::vector<Real> density;
	//	//		for (size_t p = 0; p < select_index.size(); ++p) {
	//	//			Real min_dist = 1.*10e14;
	//	//			if (std::find(sub_optima_inx.begin(), sub_optima_inx.end(), select_index[p]) == sub_optima_inx.end()) {
	//	//				//索引p所在子空间的的邻域中选中的个体之间的距离
	//	//				std::vector<size_t> compare_idx;
	//	//				auto temp_obj = pop[select_index[p]].objective();
	//	//				if (m_normalize) {
	//	//					for (size_t m = 0; m < temp_obj.size(); ++m) {
	//	//						temp_obj[m] = (temp_obj[m] - m_front_pop_range[m].first) / (m_front_pop_range[m].second - m_front_pop_range[m].first);
	//	//					}
	//	//				}
	//	//				size_t sp_idx = m_obj_tree->getRegionIdx(temp_obj);
	//	//				std::list<size_t> neighbors;
	//	//				m_obj_tree->findNeighbor(sp_idx,neighbors);
	//	//				neighbors.push_back(sp_idx);
	//	//				for (const auto ne : neighbors) {
	//	//					if (std::find(front_space_index.begin(), front_space_index.end(), ne) != front_space_index.end()) {
	//	//						for (const auto& ind : pop_att[0][0][ne]) {
	//	//							if (std::find(select_index.begin(), select_index.end(), ind) != select_index.end()) {
	//	//								compare_idx.push_back(ind);
	//	//							}
	//	//						}
	//	//					}
	//	//				}
	//	//				for (size_t k = 0; k < compare_idx.size(); ++k) {
	//	//					if (p != k) {
	//	//						auto p1 = pop[select_index[p]].objective();
	//	//						auto p2 = pop[compare_idx[k]].objective();
	//	//						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
	//	//						if (dist < min_dist) {
	//	//							min_dist = dist;
	//	//						}
	//	//					}
	//	//				}
	//	//			}
	//	//			density.push_back(min_dist);
	//	//		}
	//	//		//再去掉最拥挤的
	//	//		size_t inx = 0;
	//	//		Real dist = 1.*10e14;
	//	//		for (size_t p = 0; p < density.size(); ++p) {
	//	//			if (density[p] < dist) {
	//	//				dist = density[p];
	//	//				inx = p;
	//	//			}
	//	//		}
	//	//		select_index.erase(select_index.begin() + inx);
	//	//	}

	//	//	//当选择总数不足时,依次从前排子空间中轮流选择剩余的个体，标记选完的子空间
	//	//	std::map<size_t, std::vector<size_t>> all_ind_map= pop_att[0][0];
	//	//	
	//	//	while (select_index.size() < m_pop->size()) {
	//	//		for (const auto& sp : all_ind_map) {
	//	//			if (space_select_num[sp.first] < sp.second.size()) {//未选完的子空间
	//	//				while (select_index.size() < m_pop->size()) {
	//	//					for (const auto& ind : sp.second) {
	//	//						if (std::find(select_index.begin(), select_index.end(), ind) == select_index.end()) {
	//	//							select_index.push_back(ind);
	//	//							space_select_num[sp.first]++;
	//	//							break;
	//	//						}
	//	//					}
	//	//					break;
	//	//				}
	//	//			}
	//	//		}
	//	//	}
	//	//}
	//	//else {//前排个体数不足时，前排全选，后面的轮流选
	//	//	//先加入前排子空间个体的索引
	//	//	for (const auto &ind :pop_att[0][0]) {
	//	//		for (const auto& p : ind.second) {
	//	//			select_index.push_back(p);
	//	//		}
	//	//	}
	//	//	//再在后排子空间轮流选择
	//	//	std::vector<size_t> space_select_num(pop_att[0][1].size(),0);//每个后排子空间已选择的个数
	//	//	std::vector<std::vector<size_t>> behind_spaces;
	//	//	for (const auto& sp : pop_att[0][1]){
	//	//		behind_spaces.push_back(sp.second);
	//	//	}
	//	//	while (select_index.size() < m_pop->size()) {
	//	//		for (size_t i = 0; i < behind_spaces.size();++i) {
	//	//			if (space_select_num[i] < behind_spaces.size()) {
	//	//				size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
	//	//				if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
	//	//					select_index.push_back(behind_spaces[i][select_idx]);
	//	//					space_select_num[i]++;
	//	//				}
	//	//			}
	//	//			if (select_index.size() >= m_pop->size()) {
	//	//				break;
	//	//			}
	//	//		}
	//	//	}
	//	//}
	//	return select_index;
	//}

	//Real SPMOEA0::calSubspaceSpan(const std::vector<Real>& ideal_point, const std::vector<Real>& nadir_point, const std::vector<std::vector<Real>>& points) {
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

	//std::tuple<std::vector<size_t>,bool> SPMOEA0::selectIndSubspace(const std::vector<Real>& ref_point, const std::vector<size_t>& ind_index, int num) {
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

	//space_attach SPMOEA0::spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b) {
	//	space_attach space_all_indi;
	//	////先得到第一排个体的索引
	//	//std::vector<size_t> first_ind_index;
	//	//for (size_t i = 0; i < pop.size(); ++i) {
	//	//	if (pop[i].fitness() == 0)
	//	//		first_ind_index.emplace_back(i);
	//	//}
	//	////然后得到每个个体在哪个子空间
	//	//std::vector<size_t> space_index;
	//	//for (size_t i = 0; i < pop.size(); ++i) {
	//	//	std::vector<Real> temp;
	//	//	if (b) {
	//	//		temp = pop[i].objective();
	//	//		if (m_normalize) {
	//	//			for (size_t j = 0; j < temp.size(); ++j) {
	//	//				temp[j] = (temp[j] - m_front_pop_range[j].first) / (m_front_pop_range[j].second - m_front_pop_range[j].first);
	//	//			}
	//	//		}
	//	//	}
	//	//	else {
	//	//		temp = pop[i].variable().vect();
	//	//	}
	//	//	for (size_t j = 0; j < temp.size(); ++j) {
	//	//		if (temp[j] > 1) {
	//	//			size_t a = 1;//个体映射到子空间会越界，但是索引为最近的一个
	//	//		}
	//	//	}
	//	//	space_index.emplace_back(tree.getRegionIdx(temp));
	//	//}
	//	////得到前排子空间索引
	//	//std::vector<size_t> nondominate_space;
	//	//for (size_t i = 0; i < first_ind_index.size(); ++i) {
	//	//	size_t temp1 = space_index[first_ind_index[i]];
	//	//	if (nondominate_space.empty())
	//	//		nondominate_space.emplace_back(temp1);
	//	//	else if (std::find(nondominate_space.begin(), nondominate_space.end(), temp1) == nondominate_space.end()) {
	//	//		nondominate_space.emplace_back(temp1);
	//	//	}
	//	//}
	//	////前排子空间内的个体总数
	//	//size_t front_indi_num = 0;
	//	////得到前排子空间内包含的所有个体
	//	//std::map<size_t, std::vector<size_t>> front_indi;//key为子空间索引，value为个体索引
	//	//std::vector<size_t> redi_space_index = space_index;
	//	//for (size_t i = 0; i < nondominate_space.size(); ++i) {
	//	//	std::pair<size_t, std::vector<size_t>> temp;
	//	//	temp.first = nondominate_space[i];
	//	//	for (size_t j = 0; j < space_index.size(); ++j) {
	//	//		if (space_index[j] == temp.first) {
	//	//			temp.second.emplace_back(j);
	//	//			++front_indi_num;
	//	//		}
	//	//	}
	//	//	front_indi.insert(temp);
	//	//}
	//	//space_all_indi.emplace_back(front_indi);

	//	////非前排子空间索引
	//	//std::vector<size_t> dominated_space;
	//	////得到非前排子空间内包含的个体
	//	////std::vector<std::pair<size_t, std::vector<size_t>>> behind_indi;
	//	//std::map<size_t, std::vector<size_t>> behind_indi;
	//	//if (front_indi_num < pop.size()) {
	//	//	for (size_t i = 0; i < space_index.size(); ++i) {
	//	//		size_t temp = space_index[i];
	//	//		if (std::find(nondominate_space.begin(), nondominate_space.end(), temp) == nondominate_space.end()) {
	//	//			if (dominated_space.empty())
	//	//				dominated_space.emplace_back(temp);
	//	//			else if (std::find(dominated_space.begin(), dominated_space.end(), temp) == dominated_space.end()) {
	//	//				dominated_space.emplace_back(temp);
	//	//			}
	//	//		}
	//	//	}
	//	//	for (size_t i = 0; i < dominated_space.size(); ++i) {
	//	//		std::pair<size_t, std::vector<size_t>> temp;
	//	//		temp.first = dominated_space[i];
	//	//		for (size_t j = 0; j < space_index.size(); ++j) {
	//	//			if (space_index[j] == temp.first)
	//	//				temp.second.emplace_back(j);
	//	//		}
	//	//		behind_indi.insert(temp);
	//	//	}
	//	//}
	//	//space_all_indi.emplace_back(behind_indi);
	//	return space_all_indi;
	//}

	//pop_attach SPMOEA0::popAttach(const Population<Solution<>>& pop) {
	//	pop_attach all_indi;
	//	all_indi.emplace_back(spaceAttach(pop, getMO_HLC().getObjspaceTree(), true));
	//	all_indi.emplace_back(spaceAttach(pop, getMO_HLC().subspaceTree(), false));
	//	return all_indi;
	//}


	//void SPMOEA0::recordMetrics(Problem *pro) {
	//	/************************************/
	//	/*            性能指标计算          */
	//	/************************************/
	//	std::vector<Solution<>> temp_pop;
	//	size_t count1, count2;
	//	count1 = count2 = 0;
	//	for (size_t i = 0; i < m_pop->size(); ++i) {
	//		if (m_pop->at(i).fitness() == 0) {
	//			temp_pop.emplace_back(m_pop->at(i));
	//			count1++;
	//		}
	//		else if (m_pop->at(i).fitness() == 1) {
	//			count2++;
	//		}
	//	}
	//	m_IGD = CAST_CONOP(pro)->optima().invertGenDist(temp_pop);
	//	m_R1 = (Real)count1 / m_pop->size();
	//	m_R2 = (Real)count1 / (1 + count2) / m_pop->size();
	//	m_R3 = (Real)count1 / (1 + m_pop->size() - count1 - count2) / m_pop->size();
	//	record();//store metrics data
	//}

	void SPMOEA0::recordHisFront(Problem *pro) {
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
		if ((num_eval == 300000 && num_var < 5) || (num_eval == 500000 && num_var >= 5 && num_var < 10) || (num_eval == 700000 && num_var >= 10)) {
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

	void SPMOEA0::updatePopAge() {
		for (size_t i = 0; i < getPop().size(); i++) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				size_t temp = getPop()[i][j].timeEvaluate();
				getPop()[i][j].setTimeEvaluate(temp + 1);
			}
		}
	}

	//void SPMOEA0::generateOffspring(const std::vector<size_t>& resource, const std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>& potential_inds, Problem *pro, Random *rnd) {
	//	/************************************/
	//		/*            产生新解              */
	//		/************************************/
	//	if (!m_evolve_by_potential) {
	//		//SBX产生新解
	//		MOSP::generateOffspring(pro, rnd);
	//	}
	//	else {
	//		//再根据子空间信息确定采样，按照估计的潜力产生子代（收敛潜力和多样性潜力）
	//		//个体产生方式是否可以不用GA，按照子空间的潜力值大小分配搜索个体的数目
	//		//在每个吸引域内产生指定数量的个体：1、基于吸引域内个体生成；2、基于子空间内个体生成
	//		size_t count_offspring = 0;
	//		if (m_sample_in_basin) {//在basin内采样
	//			for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
	//				//先得到吸引域分配的解个数
	//				size_t select_num = resource[i];
	//				auto spaces = m_mo_hlc->getBasinInfo(i).m_subspace_set;
	//				//再根据根据吸引域的类型和当前解的个数产生新的解
	//				std::string basin_type = m_mo_hlc->getBasinInfo(i).flag;
	//				if (basin_type == "explore") {
	//					//在潜力最大的几个区域随机采样
	//					std::vector<size_t> temp_space = spaces;
	//					for (size_t j = 0; j < temp_space.size(); ++j) {
	//						for (size_t k = j + 1; k < temp_space.size(); ++k) {
	//							if (m_mo_hlc->getSubspaceInfo(temp_space[j]).m_potential < m_mo_hlc->getSubspaceInfo(temp_space[k]).m_potential) {
	//								auto temp = temp_space[j];
	//								temp_space[j] = temp_space[k];
	//								temp_space[k] = temp;
	//							}
	//						}
	//					}
	//					size_t count = 0;
	//					while (count < select_num) {
	//						for (size_t j = 0; j < temp_space.size(); ++j) {
	//							if (count < select_num) {
	//								//在子空间采用随机采样的方式生成一个解
	//								auto boundary_x = m_mo_hlc->subspaceTree().getBox(temp_space[j]);
	//								for (size_t k = 0; k < boundary_x.size(); ++k) {
	//									m_pop->getOffspring()[count_offspring].variable()[k] = boundary_x[k].first + (boundary_x[k].second - boundary_x[k].first) * rnd->uniform.next();
	//								}
	//								m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//								m_pop->getOffspring()[count_offspring].setCounter(0);
	//								count_offspring++;
	//								count++;
	//							}
	//							else {
	//								break;
	//							}
	//						}
	//					}
	//				}
	//				else {
	//					//使用吸引域内部当前个体产生解
	//					size_t ind_num = m_mo_hlc->getBasinInfo(i).m_current_indi.size();
	//					/*****************采用SBX产生*******************/
	//					PopSBX<> temp_pop(ind_num, pro);
	//					//初始化种群
	//					size_t count = 0;
	//					for (size_t k = 0; k < spaces.size(); ++k) {
	//						size_t num = m_mo_hlc->getSubspaceInfo(spaces[k]).m_curr_sols.size();
	//						if (num != 0) {
	//							for (size_t p = 0; p < num; ++p) {
	//								temp_pop.at(count) = *m_mo_hlc->getSubspaceInfo(spaces[k]).m_curr_sols[p];
	//								count++;
	//							}
	//						}
	//					}
	//					if (ind_num > 1) {
	//						for (size_t j = 0; j < select_num; ++j) {
	//							std::vector<size_t> p(ind_num, 0);
	//							std::vector<Real> var1;
	//							std::vector<Real> var2;
	//							Real dist = 0.;
	//							do {
	//								if (ind_num > 3) {
	//									p[0] = temp_pop.tournamentSelection(pro, rnd);
	//									p[1] = temp_pop.tournamentSelection(pro, rnd);
	//								}
	//								else {
	//									for (size_t j = 0; j < ind_num; ++j) {
	//										p[j] = j;
	//									}
	//									rnd->uniform.shuffle(p.begin(), p.end());
	//								}
	//								if (p[0] != p[1]) {
	//									var1 = temp_pop[p[0]].variable().vect();
	//									var2 = temp_pop[p[1]].variable().vect();
	//									dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
	//									if (dist < 1. * 10e-6) {
	//										size_t h_idx = std::floor(m_mo_hlc->getBasinInfo(i).m_history_inds.size() * rnd->uniform.next());
	//										temp_pop.append(*m_mo_hlc->getBasinInfo(i).m_history_inds[h_idx]);
	//										p.resize(ind_num + 1);
	//										ind_num++;
	//									}
	//								}
	//							} while (p[1] == p[0] || dist <= 1. * 10e-6);
	//							//判断交叉之后子代时候与父代相同，相同则继续
	//							do {
	//								temp_pop.crossover(p[0], p[1], m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
	//								temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
	//								//子代越界处理
	//								//for (size_t j = 0; j < m_pop_var_range.size(); ++j) {
	//								//	auto range = m_pop_var_range[j];
	//								//	auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
	//								//	while (var1 < range.first) {
	//								//		m_pop->getOffspring()[count_offspring].variable().vect()[j] = 2 * range.first - var1;
	//								//		var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
	//								//	}
	//								//	while (var1 > range.second) {
	//								//		m_pop->getOffspring()[count_offspring].variable().vect()[j] = 2 * range.second - var1;
	//								//		var1 = m_pop->getOffspring()[count_offspring].variable().vect()[j];
	//								//	}
	//								//	//判断是否越界
	//								//	auto tmp = m_pop->getOffspring()[count_offspring].variable().vect()[j];
	//								//	int isnan(tmp);
	//								//	if (!isnormal(tmp)) {
	//								//		size_t a = 1;
	//								//	}
	//								//	if (tmp<range.first || tmp>range.second) {
	//								//		size_t a = 1;
	//								//	}
	//								//}
	//								repairSol(m_pop->getOffspring()[count_offspring].variable().vect(), pro, rnd);
	//							} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[0]].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[1]].variable().vect());

	//							m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//							m_pop->getOffspring()[count_offspring].setCounter(0);
	//							count_offspring++;
	//						}
	//					}
	//					else {
	//						//采用局部变异的方式，先找出个体所在子空间
	//						size_t idx = 0;
	//						for (size_t j = 0; j < spaces.size(); ++j) {
	//							if (!m_mo_hlc->getSubspaceInfo(spaces[j]).m_curr_sols.empty()) {
	//								idx = j;
	//							}
	//						}
	//						auto temp_x = m_mo_hlc->getSubspaceInfo(spaces[idx]).m_curr_sols.back()->variable().vect();
	//						auto boundary_x = m_mo_hlc->subspaceTree().getBox(spaces[idx]);
	//						for (size_t j = 0; j < select_num; ++j) {
	//							auto temp_off = localMutationSol(temp_x, boundary_x, pro, rnd);
	//							m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
	//							m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//							m_pop->getOffspring()[count_offspring].setCounter(0);
	//							count_offspring++;
	//						}
	//					}
	//				}
	//			}
	//		}
	//		else {//在子空间采样
	//			for (size_t i = 0; i < m_mo_hlc->numBasin(); ++i) {
	//				//在吸引域中，根据子空间的搜索潜力决定在哪个子空间采样
	//				std::string basin_type = m_mo_hlc->getBasinInfo(i).flag;
	//				size_t select_num = resource[i];
	//				std::vector<size_t> space_idx = m_mo_hlc->getBasinInfo(i).m_subspace_set;
	//				std::vector<Real> subspace_potential;
	//				Real sum = 0.;
	//				for (size_t j = 0; j < space_idx.size(); ++j) {
	//					subspace_potential.push_back(m_mo_hlc->getSubspaceInfo(space_idx[j]).m_potential);
	//					sum += subspace_potential.back();
	//				}
	//				for (auto& i : subspace_potential) {
	//					i = i / sum;
	//				}
	//				auto select_probability = subspace_potential;
	//				for (size_t i = 0; i < subspace_potential.size(); ++i) {
	//					Real sum = 0.;
	//					for (size_t j = 0; j <= i; ++j) {
	//						sum += subspace_potential[j];
	//					}
	//					select_probability[i] = sum;
	//				}
	//				for (size_t i = 0; i < select_num; ++i) {
	//					Real rand_num = rnd->uniform.next();
	//					size_t idx;
	//					for (size_t j = 0; j < select_probability.size(); ++j) {
	//						if (rand_num <= select_probability[j]) {
	//							idx = j;
	//							break;
	//						}
	//					}
	//					//选择采样的子空间为space_idx[idx],并在里面进行采样
	//					auto boundary_x = m_mo_hlc->subspaceTree().getBox(space_idx[idx]);
	//					if (basin_type == "explore" || m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.empty()) {//随机采样的方式
	//						for (size_t k = 0; k < boundary_x.size(); ++k) {
	//							m_pop->getOffspring()[count_offspring].variable()[k] = boundary_x[k].first + (boundary_x[k].second - boundary_x[k].first) * rnd->uniform.next();
	//						}
	//					}
	//					else {//基于历史解生成，有问题，生成解
	//						if (m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.size() > 1) {
	//							size_t pop_num = m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.size();
	//							/*****************采用SBX产生*******************/
	//							PopSBX<> temp_pop(pop_num, pro);
	//							//初始化种群
	//							for (size_t j = 0; j < pop_num; ++j) {
	//								temp_pop.at(j) = *m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols[j];
	//							}
	//							std::vector<size_t> p(pop_num, 0);
	//							std::vector<Real> var1;
	//							std::vector<Real> var2;
	//							Real dist = 0.;
	//							do {
	//								if (pop_num > 3) {
	//									p[0] = temp_pop.tournamentSelection(pro, rnd);
	//									p[1] = temp_pop.tournamentSelection(pro, rnd);
	//								}
	//								else {
	//									for (size_t j = 0; j < pop_num; ++j) {
	//										p[j] = j;
	//									}
	//									rnd->uniform.shuffle(p.begin(), p.end());
	//								}
	//								if (p[0] != p[1]) {
	//									var1 = temp_pop[p[0]].variable().vect();
	//									var2 = temp_pop[p[1]].variable().vect();
	//									dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
	//									if (dist < 1. * 10e-6) {
	//										size_t h_idx = std::floor(m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds.size() * rnd->uniform.next());
	//										temp_pop.append(*m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds[h_idx]);
	//										p.resize(pop_num + 1);
	//										pop_num++;
	//									}
	//								}
	//							} while (p[1] == p[0] || dist <= 1. * 10e-6);
	//							do {
	//								temp_pop.crossover(p[0], p[1], m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
	//								temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
	//								//子代越界处理
	//								/*for (size_t k = 0; k < m_pop_var_range.size(); ++k) {
	//									auto range = m_pop_var_range[k];
	//									auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//									while (var1 < range.first) {
	//										m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.first - var1;
	//										var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//									}
	//									while (var1 > range.second) {
	//										m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.second - var1;
	//										auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//									}
	//								}*/
	//								repairSol(m_pop->getOffspring()[count_offspring].variable().vect(), pro, rnd);
	//							} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[0]].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[p[1]].variable().vect());

	//							/*****************采用DE产生*******************/
	//							//PopDE<> temp_pop(pop_num,pro);
	//							////初始化种群
	//							////交叉变异

	//						}
	//						else {//当只有一个解时，在解的邻域产生,每一维的3 sigma
	//							auto temp_x = m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_curr_sols.back()->variable().vect();
	//							auto boundary_x = m_mo_hlc->subspaceTree().getBox(space_idx[idx]);
	//							auto temp_off = localMutationSol(temp_x, boundary_x, pro, rnd);
	//							m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
	//						}
	//					}
	//					m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//					m_pop->getOffspring()[count_offspring].setCounter(0);
	//					count_offspring++;
	//				}
	//			}
	//		}
	//		//在稀疏点处采样
	//		size_t num_basin = m_mo_hlc->numBasin();
	//		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
	//		std::vector<std::vector<size_t>> poten_inds = potential_inds.second;
	//		for (size_t i = num_basin; i < resource.size(); ++i) {
	//			if (i >= num_basin && i < num_basin + num_obj) {
	//				//依次得到每个点所在的子空间，然后向延伸方向采样
	//				auto pop_idx = poten_inds[i - num_basin][0];
	//				auto var_position = m_pop->at(pop_idx).variable().vect();
	//				auto obj_position = m_pop->at(pop_idx).objective();
	//				auto var_index = m_mo_hlc->subspaceTree().getRegionIdx(var_position);
	//				auto boundary_x = m_mo_hlc->subspaceTree().getBox(var_index);
	//				auto obj_index = m_mo_hlc->getObjspaceTree().getRegionIdx(obj_position);
	//				auto ind_index = getObjRegionInfo(obj_index).ind_idx;
	//				if (ind_index.size() > 1) {
	//					std::vector<Real> obj1 = m_pop->at(pop_idx).objective();;
	//					std::vector<Real> obj2 = m_pop->at(ind_index[0]).objective();;
	//					Real dist = 0.;
	//					size_t count = 0;
	//					std::vector<Real> ind_var;
	//					do {
	//						rnd->uniform.shuffle(ind_index.begin(), ind_index.end());
	//						obj2 = m_pop->at(ind_index[0]).objective();
	//						dist = euclideanDistance(obj1.begin(), obj1.end(), obj2.begin());
	//						if (dist < 1. * 10e-6) {
	//							count++;
	//						}
	//						if (count > 20) {
	//							break;
	//						}
	//					} while (ind_index[0] == pop_idx || dist < 1. * 10e-7);
	//					ind_var = m_pop->at(ind_index[0]).variable().vect();
	//					if (count > 20) {
	//						size_t his_sol_num = m_mo_hlc->getSubspaceInfo(var_index).m_history_inds.size();
	//						std::vector<Real> var1;
	//						do {
	//							size_t select_idx = std::floor(his_sol_num * rnd->uniform.next());
	//							var1 = (*m_mo_hlc->getSubspaceInfo(var_index).m_history_inds[select_idx]).variable().vect();
	//							dist = euclideanDistance(var1.begin(), var1.end(), var_position.begin());
	//						} while (dist < 1. * 10e-4);
	//						ind_var = var1;
	//					}
	//					auto temp_off = vectorMutationSol(var_position, ind_var, pro, rnd);
	//					m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
	//					m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//					m_pop->getOffspring()[count_offspring].setCounter(0);
	//					count_offspring++;
	//				}
	//				//if (m_mo_hlc->getSubspaceInfo(var_index).m_curr_sols.size() > 1) {
	//				//	//方向采样
	//				//	auto ind_idx = m_mo_hlc->getSubspaceInfo(var_index).m_curr_sols_idx;
	//				//	for (size_t k = 0; k < resource[i]; ++k) {
	//				//		//任意找一个非边界点,但不和边界点重合
	//				//		std::vector<Real> var1;
	//				//		std::vector<Real> var2;
	//				//		do {
	//				//			rnd->uniform.shuffle(ind_idx.begin(), ind_idx.end());
	//				//			var1 = temp_pop[p[0]].variable().vect();
	//				//			var2 = temp_pop[p[1]].variable().vect();
	//				//			dist = euclideanDistance(var1.begin(), var1.end(), var2.begin());
	//				//			if (dist < 1. * 10e-6) {
	//				//				size_t h_idx = std::floor(m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds.size() * rnd->uniform.next());
	//				//				temp_pop.append(*m_mo_hlc->getSubspaceInfo(space_idx[idx]).m_history_inds[h_idx]);
	//				//				p.resize(pop_num + 1);
	//				//				pop_num++;
	//				//			}
	//				//		} while (ind_idx[0] == pop_idx);
	//				//		//向量计算，步长设置
	//				//		auto select_ind = m_pop->at(ind_idx[0]).variable().vect();
	//				//		auto temp_off = vectorMutationSol(var_position,select_ind,rnd);
	//				//		m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
	//				//		m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//				//		m_pop->getOffspring()[count_offspring].setCounter(0);
	//				//		count_offspring++;
	//				//	}
	//				//}
	//				else {
	//					//局部变异
	//					for (size_t k = 0; k < resource[i]; ++k) {
	//						auto temp_off = localMutationSol(var_position, boundary_x, pro, rnd);
	//						m_pop->getOffspring()[count_offspring].variable().vect() = temp_off;
	//						m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//						m_pop->getOffspring()[count_offspring].setCounter(0);
	//						count_offspring++;
	//					}
	//				}
	//			}
	//			else {
	//				//先得到两个个体及其所在的子空间
	//				//auto ind1 = m_pop->at((*poten_inds[i-num_basin+num_obj])[0]).variable().vect();
	//				//auto ind2 = m_pop->at((*poten_inds[i - num_basin + num_obj])[1]).variable().vect();
	//				/*****************采用SBX产生*******************/
	//				PopSBX<> temp_pop(2, pro);
	//				temp_pop.setEta(2, 2);
	//				//初始化种群
	//				for (size_t j = 0; j < 2; ++j) {
	//					temp_pop.at(j) = m_pop->at(poten_inds[i - num_basin][j]);
	//				}
	//				for (size_t j = 0; j < resource[i]; ++j) {
	//					do {
	//						temp_pop.crossover(0, 1, m_pop->getOffspring()[count_offspring], m_pop->getParent()[count_offspring], pro, rnd);
	//						temp_pop.mutate(m_pop->getOffspring()[count_offspring], pro, rnd);
	//						//子代越界处理
	//						for (size_t k = 0; k < m_pop_var_range.size(); ++k) {
	//							auto range = m_pop_var_range[k];
	//							auto var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//							while (var1 < range.first) {
	//								m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.first - var1;
	//								var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//							}
	//							while (var1 > range.second) {
	//								m_pop->getOffspring()[count_offspring].variable().vect()[k] = 2 * range.second - var1;
	//								var1 = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//							}
	//							////越界判断
	//							//auto tmp = m_pop->getOffspring()[count_offspring].variable().vect()[k];
	//							//if (!isnormal(tmp)) {
	//							//	size_t a = 1;
	//							//}
	//							//if (tmp<m_pop_var_range[k].first || tmp>m_pop_var_range[k].second) {
	//							//	size_t a = 1;
	//							//}
	//						}
	//						repairSol(m_pop->getOffspring()[count_offspring].variable().vect(), pro, rnd);
	//					} while (m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[0].variable().vect() || m_pop->getOffspring()[count_offspring].variable().vect() == temp_pop[1].variable().vect());
	//					m_pop->getOffspring()[count_offspring].setId(count_offspring);
	//					m_pop->getOffspring()[count_offspring].setCounter(0);
	//					count_offspring++;
	//				}
	//			}
	//		}
	//	}
	//}
}
