#include "SPMOEA1_0.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA1_0::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_0::run_() {
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

	void SPMOEA1_0::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA1_0::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(getPop().size() + 1);//���һ��Ϊ��ʷ��
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

	void SPMOEA1_0::initiVarSpace(Problem *pro) {
		/*************************************/
		/*      MO_HLC�����ռ仮�ֳ�ʼ��     */
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

	void SPMOEA1_0::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
		//�ڸ����ӿռ�����һ������Ⱥ
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
		Population<Solution<>> new_pop;
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
			SPMOEA::NDSort(temp_pop);
			//��ʼ���Ӵ�
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring()[j] = temp_pop[j];
				temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
				new_pop.append(temp_pop[j]);
			}
			/*temp_pop.setRate(getCr(), getMr());
			temp_pop.setEta(getCeta(), getMeta());*/
			temp_pop.setPopState("active");
			//��ʼ������Ⱥ����ʷ��ÿ��ǰ�ؽ�
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
			getMO_HLC().getSubspaceInfo(i).m_feature_state = 1;
			getMO_HLC().getSubspaceInfo(i).idx_cluster.push_back(i);
		}
		//������ʷ��Ϣ
		SPMOEA::NDSort(new_pop);
		updateHistoryInfo(new_pop, pro);
		SPMOEA::updateObjRange(new_pop, pro);
		SPMOEA::updateNewPop(new_pop);
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();
		//�����ӿռ仮�־���
		clusterSubspace();
		//�����ӿռ���Ϣ
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				Population<Solution<>> temp_off_pop;
				for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); ++j) {
					temp_off_pop.append(getPop()[i].getOffspring()[j]);
				}
				getMO_HLC().updateSubspaceInfo(temp_off_pop, pro, rnd);
			}
		}
		updateVarSpace(pro, rnd);
		updateEE(space_num * pop_num, 0);
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
	}

	int SPMOEA1_0::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		/********************************************************************************
							   ���ÿ������Ⱥ״̬���Ӷ��ж��ӿռ�ṹ
		********************************************************************************/
		std::vector<Real> front_in_bound_ratio;
		calFrontBoundRatio(front_in_bound_ratio);
		/********************************************************************************
							   ��������Ⱥǰ���ӿռ���������������Դ
		********************************************************************************/
		//����ʹ������Ⱥ��Ϣ���ڲ���ȫ��Ϣ�Ƚϣ��ɷ�ʹ���ӿռ�ȫ������Ч��Ϣ�����ӿռ�ıȽϣ�
		//����������Ⱥ�����������
		//����ǰ���front_in_bound_ratioȷ��ÿ������Ⱥ��������Դ��̽��������
		std::vector<Real> pop_exploit_ratio;
		std::vector<size_t> assign_pop_resource;
		PopResourceAssign(pop_exploit_ratio, assign_pop_resource, front_in_bound_ratio, pro);
		/********************************************************************************
											 ����Ⱥ�ݻ�
		********************************************************************************/
		//1.�Ӵ���̽���뿪���ı�����2.�Ӵ���̽���뿪���ķ�ʽ
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			bool m_search_balance = false;//ѡ������Ⱥ�Ƿ�E&Eƽ��
			if (!m_search_balance) {
				for (auto& i : pop_exploit_ratio) {
					i = 1.;
				}
			}
			bool m_assign_resource = true;//����Ⱥ�Ƿ������Դ����
			if (!m_assign_resource) {
				for (auto& i : assign_pop_resource) {
					i = getPopsize();
				}
			}
			generateOffspring(pro, rnd, assign_pop_resource, pop_exploit_ratio, 0);
		}
		else {
			//generating solutions by learning methods
			//multi-popution and prediction
			//��������Ԥ�⣬����Ԥ�����
			//����Ⱥ�ݻ���ʷ�켣
		}
		Population<Solution<>> offspring_pop;
		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
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
			}
		}
		/**********************************************************************************
					����Ƿ�ϸ���ӿռ�:���������Ⱥ���������ռ��ǰ���ӿռ�ϸ��
		***********************************************************************************/
		spaceSubdivision(pro, rnd);
		clusterSubspace();
		/**********************************************************************************
							   �ۺϸ�������Ⱥ��Ϣ�������ӿռ���Ϣ
		**********************************************************************************/
		SPMOEA::NDSort(offspring_pop);
		//updateVarSpace(pro, rnd);//ʹ������Ⱥ�Ӵ������ӿռ���Ϣ
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//ʹ����ʷ���з�֧������archive
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
										 ��Ⱥ���£�����Ⱥ�ڲ�������̭ѡ��
		*********************************************************************************************/
		for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
				}
				auto& bound = getPop()[i].getPopSearchBox();
				getPop()[i].envirSelection(pro, rnd, 0);
				//getPop()[i].iteration()++;

				updatePopDistribute(i, pro);
				updatePopdist(i);
			}
		}
		/*********************************************************************************************
												 ��¼������Ϣ
		**********************************************************************************************/
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		return tag;
	}

	void SPMOEA1_0::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			auto pop_state = getPop()[i].getPopState();
			auto search_box = getPop()[i].getPopSearchBox();
			if (pop_state == "active") {
				//��Ⱥ��ʷǰ�ؽ���������ĸ���
				size_t in_count = 0;
				auto& pop_his_front_sols = getPop()[i].getPopHisFrontSols();
				for (size_t j = 0; j < pop_his_front_sols.size(); ++j) {
					if (indInRange(search_box, pop_his_front_sols[j]->variable().vect())) {
						in_count++;
					}
				}
				Real in_ratio = (Real)in_count / pop_his_front_sols.size();
				front_bound_ratio.push_back(in_ratio);
				if (in_ratio < 0.1) {//����
					getPop()[i].setPopState("sleep");
				}
			}
			else {
				front_bound_ratio.push_back(front_bound_ratio.back());
			}
		}
	}

	bool SPMOEA1_0::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA1_0::updateSubPopSpace(size_t pop_inx, Problem *pro, Random *rnd) {
		//��������Ⱥ����Ϣ����������Ⱥ������������ӿռ�����ֵ
		//�����ӿռ�������Ⱥ���������ڵ����rankֵ
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
			//������Ⱥ�������ӿռ��������ֵ
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
			//����û�и�����ӿռ䣬��������õ�����ֵ
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

	void SPMOEA1_0::updateVarSpace(Problem *pro, Random *rnd) {
		/*for (size_t i = 0; i < getPop().size(); ++i) {
			if (getPop()[i].getPopState() == "active") {
				Population<Solution<>> temp_off_pop;
				for (size_t j = 0; j < getPop()[i].getOffspring().size()- getPop()[i].size(); ++j) {
					temp_off_pop.append(getPop()[i].getOffspring()[j]);
				}
				getMO_HLC().updateSubspaceInfo(temp_off_pop, pro, rnd);
			}
		}*/

		////�����ӿռ�������Ⱥ���������ڵ����rankֵ
		//for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
		//	if (getPop()[i].getPopState() == "active") {
		//		std::vector<std::shared_ptr<Solution<>>> temp_pop1;
		//		for (size_t j = 0; j < getMO_HLC().getCluster(i)[0].size(); ++j) {
		//			for (size_t k = 0; k < getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_represent_sol.size(); ++k) {
		//				temp_pop1.emplace_back(getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_represent_sol[k]);
		//			}
		//		}
		//		std::vector<std::vector<Real>*> objs;
		//		for (size_t j = 0; j < temp_pop1.size(); ++j) {
		//			objs.emplace_back(&temp_pop1[j]->objective());
		//		}
		//		std::vector<int> rank;
		//		int num_layer = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		//		//getMO_HLC().setLayerNum(num_layer);
		//		for (size_t j = 0; j < temp_pop1.size(); ++j) {
		//			temp_pop1[j]->setFitness(rank[j]);
		//		}
		//		//�����ӿռ��������ֵ
		//		for (size_t j = 0; j < getMO_HLC().getCluster(i)[0].size(); ++j) {
		//			getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_sub_best_rank = INT16_MAX;
		//		}
		//		for (size_t j = 0; j < temp_pop1.size(); ++j) {
		//			auto temp_sol = temp_pop1[j]->variable().vect();
		//			auto inx = getMO_HLC().subspaceTree().getRegionIdx(temp_sol);
		//			if (temp_pop1[j]->fitness() < getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank) {
		//				getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank = temp_pop1[j]->fitness();
		//			}
		//		}
		//		//����û�и�����ӿռ䣬��������õ�����ֵ
		//		std::vector<size_t> cluster_spaces = getMO_HLC().getCluster(i)[0];
		//		for (size_t j = 0; j < cluster_spaces.size(); ++j) {
		//			if (getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_best_rank == INT16_MAX) {
		//				auto neigh_space = getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_neighbors;
		//				std::vector<int> neigh_ranks;
		//				/*std::vector<size_t> neigh_space;
		//				for (auto iter = temp_neigh_space.begin(); iter != temp_neigh_space.end(); ++iter) {
		//					neigh_space.push_back(*iter);
		//				}*/
		//				for (auto iter = neigh_space.begin(); iter != neigh_space.end(); ++iter) {
		//					size_t inx = *iter;
		//					if (std::find(cluster_spaces.begin(), cluster_spaces.end(), inx) != cluster_spaces.end()) {
		//						if (getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank < INT16_MAX) {
		//							neigh_ranks.push_back(getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank);
		//						}
		//					}
		//				}
		//				/*for (size_t k = 0; k < neigh_space.size(); ++k) {
		//					if (std::find(getMO_HLC().getCluster(i)[0].begin(), getMO_HLC().getCluster(i)[0].end(), neigh_space[k]) != getMO_HLC().getCluster(i)[0].end()) {
		//						if (getMO_HLC().getSubspaceInfo(neigh_space[k]).m_sub_best_rank < INT16_MAX) {
		//							neigh_ranks.push_back(getMO_HLC().getSubspaceInfo(neigh_space[k]).m_sub_best_rank);
		//						}
		//					}
		//				}*/
		//				if (neigh_ranks.empty()) {
		//					getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_sub_best_rank = num_layer;
		//				}
		//				else if (neigh_ranks.size() == 1) {
		//					getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_sub_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
		//				}
		//				else {
		//					int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
		//					int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
		//					getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(i)[0][j]).m_sub_best_rank = std::ceil((min_rank + max_rank)/2);
		//				}
		//			}
		//		}
		//	}
		//}

		//�����ӿռ���������ֵ��ȡ�ӿռ�Ĵ�����廹��ǰ�ظ��壿
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
		//�����ӿռ��������ֵ
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
		//����û�и�����ӿռ䣬��������õ�����ֵ
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

	void SPMOEA1_0::PopResourceAssign(std::vector<Real>& pop_exploit_ratio, std::vector<size_t>& assign_pop_resource, std::vector<Real>& front_bound_ratio, Problem *pro) {
		//������Ⱥ����ʷǰ�ؽ��ڱ߽��ڵı������������Դ
		size_t max_pop_size = getPopsize();
		for (size_t i = 0; i < front_bound_ratio.size(); ++i) {
			assign_pop_resource.push_back(std::ceil(front_bound_ratio[i] * max_pop_size));
			pop_exploit_ratio.push_back(1.);
		}
		//std::vector<Real> cluster_front_volume(getPop().size(),0.);//����Ⱥǰ���ӿռ����
		//for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//	auto id_clu = getMO_HLC().getSubspaceInfo(i).idx_cluster.front();
		//	if (getPop()[id_clu].getPopState() == "active") {
		//		if (getMO_HLC().getSubspaceInfo(i).m_sub_best_rank == 0) {
		//			auto v = getMO_HLC().subspaceTree().getBoxVolume(i);
		//			cluster_front_volume[id_clu] += v;
		//		}
		//	}
		//	else {
		//		cluster_front_volume[id_clu] += 0.;
		//	}
		//}
		//auto max_v = *std::max_element(cluster_front_volume.begin(),cluster_front_volume.end());
		/*for (size_t i = 0; i < cluster_front_volume.size(); ++i) {
			if (cluster_front_volume[i] == 0.) {
				assign_pop_resource.push_back(0);
			}
			else {
				assign_pop_resource.push_back(std::ceil(cluster_front_volume[i]/max_v*max_pop_size));
			}
			pop_exploit_ratio.push_back(1.);
		}*/
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA1_0::generateOffspring(Problem *pro, Random *rnd, const std::vector<size_t>& pop_resource, const std::vector<Real>& pop_exploit_ratio, int type) {
		//size_t num_explore = 0;
		//size_t num_exploit = 0;
		//for (size_t i = 0; i < getPop().size(); ++i) {
		//	size_t explore_num = std::ceill(pop_resource[i] * (1 - pop_exploit_ratio[i]));
		//	size_t exploit_num = pop_resource[i] - explore_num;
		//	//if (exploit_num == 0) {//ʹ����Ⱥ������һ��������̽����
		//	//	exploit_num++;
		//	//	explore_num--;
		//	//}
		//	//if (explore_num == 0) {//ʹ����Ⱥ������һ��������̽����
		//	//	explore_num++;
		//	//	exploit_num--;
		//	//}
		//	num_explore += explore_num;
		//	num_exploit += exploit_num;
		//	//�Ӵ�������������
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
		//	//������Ⱥǰ�Ÿ���俪��,������ʽѡ�����������֮�佻��
		//	auto search_bound = CAST_CONOP(pro)->boundary();
		//	for (size_t j = 0; j < exploit_num; ++j) {
		//		std::vector<size_t> p(2);
		//		do {//���ظ�����
		//			do {
		//				p[0] = getPop()[i].tournamentSelection(pro, rnd);
		//				p[1] = getPop()[i].tournamentSelection(pro, rnd);
		//			} while (p[1] == p[0]);
		//			getPop()[i].crossover(p[0], p[1], getPop()[i].getOffspring()[j], getPop()[i].getOffspring()[j + 1], pro, rnd);
		//			getPop()[i].mutate(getPop()[i].getOffspring()[j], pro, rnd);
		//			getPop()[i].mutate(getPop()[i].getOffspring()[j + 1], pro, rnd);
		//			//�Ӵ�Խ�紦��
		//			repairSol(getPop()[i].getOffspring()[j].variable().vect(), search_bound, rnd);
		//			repairSol(getPop()[i].getOffspring()[j + 1].variable().vect(), search_bound, rnd);
		//		} while (getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//			getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[1]].variable().vect() || \
		//			getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//			getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[1]].variable().vect());

		//		getPop()[i].getOffspring()[j].setCounter(0);
		//		getPop()[i].getOffspring()[j + 1].setCounter(0);

		//		//ÿ��ȷ��һ���⣬��2���Ӵ���ѡ������Ÿ�����ģ�����������廥��֧�䣬ѡ���븸����С��������
		//		Dominance p_ship = objectiveCompare(getPop()[i][p[0]].objective(), getPop()[i][p[1]].objective(), CAST_CONOP(pro)->optimizeMode());
		//		auto point1 = getPop()[i].getOffspring()[j].variable().vect();
		//		auto point2 = getPop()[i].getOffspring()[j + 1].variable().vect();
		//		if (p_ship == Dominance::kDominant) {
		//			//ѡ��p[0]����
		//			Real s1 = euclideanDistance(point1.begin(), point1.end(), getPop()[i][p[0]].variable().vect().begin());
		//			Real s2 = euclideanDistance(point2.begin(), point2.end(), getPop()[i][p[0]].variable().vect().begin());
		//			if (s1 > s2) {
		//				getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//			}
		//		}
		//		else if (p_ship == Dominance::kDominated) {
		//			//ѡ��p[1]����
		//			Real s1 = euclideanDistance(point1.begin(), point1.end(), getPop()[i][p[1]].variable().vect().begin());
		//			Real s2 = euclideanDistance(point2.begin(), point2.end(), getPop()[i][p[1]].variable().vect().begin());
		//			if (s1 > s2) {
		//				getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//			}
		//		}
		//		else if (p_ship == Dominance::kNonDominated) {
		//			//ѡ�븸������С�����С��
		//			Real s11 = euclideanDistance(point1.begin(), point1.end(), getPop()[i][p[0]].variable().vect().begin());
		//			Real s12 = euclideanDistance(point1.begin(), point1.end(), getPop()[i][p[1]].variable().vect().begin());
		//			Real s1 = s11 < s12 ? s11 : s12;
		//			Real s21 = euclideanDistance(point2.begin(), point2.end(), getPop()[i][p[0]].variable().vect().begin());
		//			Real s22 = euclideanDistance(point2.begin(), point2.end(), getPop()[i][p[1]].variable().vect().begin());
		//			Real s2 = s21 < s22 ? s21 : s22;
		//			if (s1 > s2) {
		//				getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//			}
		//		}
		//		/*Real rnd_num = rnd->uniform.next();
		//		if (rnd_num > 0.5) {
		//			getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//		}*/
		//	}

		//	//������Ⱥ�����������ڵ��ӿռ��ڽ���̽��
		//	//����̽�������Ǽ����ӿռ�ĸ�����
		//	bool random_sample = true;
		//	std::vector<size_t> sample_spaces;
		//	if (random_sample) {//�������
		//		for (size_t j = 0; j < explore_num; ++j) {
		//			size_t ind_inx = std::floor(total_spaces.size() * rnd->uniform.next());
		//			sample_spaces.push_back(total_spaces[ind_inx]);
		//		}
		//	}
		//	else {//���ڸ����ʲ���
		//		//ȷ��ÿ���ӿռ�Ļ��ַ����������ӿռ�Ĵ�С�����Ժ��ӿռ�����ʷ��ĸ���
		//		std::vector<Real> space_density;
		//		for (size_t j = 0; j < total_spaces.size(); ++j) {
		//			auto volume = getMO_HLC().subspaceTree().getBoxVolume(total_spaces[j]);
		//			space_density.push_back(getMO_HLC().getSubspaceInfo(total_spaces[j]).m_history_inds.size() / volume);
		//		}
		//		//���ܶ���С���ӿռ����
		//		size_t select_count = 0;
		//		while (select_count < explore_num) {
		//			//����͵ĸ����ʵ��ӿռ�
		//			size_t min_inx = std::distance(space_density.begin(), std::min_element(space_density.begin(), space_density.end()));
		//			sample_spaces.push_back(total_spaces[min_inx]);
		//			//�����ӿռ�ĸ�����
		//			auto volume = getMO_HLC().subspaceTree().getBoxVolume(total_spaces[min_inx]);
		//			space_density[min_inx] = space_density[min_inx] + 1. / volume;
		//			select_count++;
		//		}
		//	}
		//	//���ӿռ��������������ǻ������и�������½�
		//	for (size_t j = 0; j < sample_spaces.size(); ++j) {
		//		auto bound = getMO_HLC().subspaceTree().getBox(sample_spaces[j]);
		//		for (size_t k = 0; k < bound.size(); ++k) {
		//			getPop()[i].getOffspring()[exploit_num + j].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
		//		}
		//	}
		//}
		//updateEE(num_explore, num_exploit);
	}

	void SPMOEA1_0::spaceSubdivision(Problem *pro, Random *rnd) {
		//��������Ⱥ�����������ڵ�����Ⱥǰ���ӿռ�
		Real total_volume = getVarSpaceVolume();
		size_t num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		for (size_t i = 0; i < getPop().size(); ++i) {
			auto bound = getPop()[i].getSearchRange();
			if (getPop()[i].getPopState() == "active") {
				//��ʹ���Ӵ��������ӿռ������ֵ�����ҳ���ǰ���ӿռ�
				updateSubPopSpace(i, pro, rnd);
				std::vector<size_t> front_space;
				auto cluster_spaces = getMO_HLC().getCluster(i)[0];
				for (size_t j = 0; j < cluster_spaces.size(); ++j) {
					if (getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_best_rank == 0) {
						front_space.push_back(cluster_spaces[j]);
					}
				}
				//ϸ���ӿռ�
				for (auto& sp : front_space) {
					auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
					Real min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
					min_num = 20;
					if (space_volume / total_volume > std::pow(1. / min_num, num_var)) {
						//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
						int dim = findSplitDim(sp, pro);
						auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
						Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
						splitSpace(sp, 2, dim, pos, false, pro, rnd);
					}
				}
			}
		}

		//�ӿռ����ռ��
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
	}

	void SPMOEA1_0::clusterSubspace() {
		//ÿ�ξ���ֻ����һ�㣬����������ص�
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		std::map<int, std::vector<size_t>> cluster;
		for (size_t i = 0; i < space_num; ++i) {
			int cluster_id = getMO_HLC().getSubspaceInfo(i).idx_cluster.front();
			if (cluster[cluster_id].empty()) {
				std::vector<size_t> temp;
				temp.push_back(i);
				cluster.insert(std::make_pair(cluster_id, temp));
			}
			cluster[cluster_id].push_back(i);
		}
		for (auto& clu : cluster) {
			std::vector<std::vector<size_t>> temp_cluster;
			temp_cluster.emplace_back(clu.second);
			getMO_HLC().getClusters().emplace_back(temp_cluster);
		}
	}

	bool SPMOEA1_0::subspaceLink(size_t inx1, size_t inx2) {
		bool flag = false;
		//�����ӿռ��е�ǰ�ؽ�ķֲ�
		auto ind1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
		auto ind2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
		////���ӿռ�߽��ƽ������
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
		//����������е���̾���ֵ
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < ind1.size(); ++i) {
			for (size_t j = 0; j < ind2.size(); ++j) {
				auto dist = euclideanDistance(ind1[i]->variable().vect().begin(), ind1[i]->variable().vect().end(), ind2[j]->variable().vect().begin());
				min_v = min_v < dist ? min_v : dist;
			}
		}
		//������и��Եĵ�����ֵ
		Real max_dist1 = 0.;
		Real max_dist2 = 0;
		std::vector<Real> dim_span1;
		std::vector<Real> dim_span2;
		//Ѱ��ÿһά������
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

	void SPMOEA1_0::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA1_0::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem *pro, Random *rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateVarSpaceInfo(temp_pop, pro, rnd);
	}

	void SPMOEA1_0::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA1_0::findSplitDim(int inx, Problem *pro) {
		//��������һά��ƽ�����ߵ��ӿռ�ǰ�ؽ�Ĳ�������ȷ�����ֵ�ά��
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

		//ʹ�ÿ��ȷ���ָ��ά��
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		return dim;
	}

	void SPMOEA1_0::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA1_0::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            ����ָ�����          */
		/************************************/
		//ʹ��archive��������ָ��
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

	void SPMOEA1_0::initiObjSpace(Problem *pro) {

	}
}