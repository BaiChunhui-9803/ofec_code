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
			m_solution.resize(getPop().size() + 1);//���һ��Ϊ��ʷ��
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

	void SPMOEA1_0_1::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//�����ӿռ�Ԥ���������ݱȽϽ�����ɶ���Ⱥ
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
		//��ʼ���Ӵ�
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
		}
		/*temp_pop.setRate(getCr(), getMr());
		temp_pop.setEta(getCeta(), getMeta());*/
		getPop().append(temp_pop);

		//�����ӿռ���Ϣ
		updateVarSpaceInfo(temp_pop, pro, rnd);
		SPMOEA::updateHistoryInfo(temp_pop, pro);
		SPMOEA::updateArchive(archiveNum(), pro);
		SPMOEA::updateNewPop(temp_pop);
		SPMOEA::updateObjRange(temp_pop, pro);
		updateObjSpace();
		//�����ӿռ�����ֵ����
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
				 ���ݸ�����Ⱥ�ıȽϣ����ø�����Ⱥ�Ӵ���̽�������������Ӵ��ĸ�����
		********************************************************************************/
		//����ʹ������Ⱥ��Ϣ���ڲ���ȫ��Ϣ�Ƚϣ��ɷ�ʹ���ӿռ�ȫ������Ч��Ϣ�����ӿռ�ıȽϣ�
		//����������Ⱥ�����������
		std::vector<Real> pop_exploit_ratio;
		std::vector<size_t> assign_pop_resource;
		PopResourceAssign(pop_exploit_ratio, assign_pop_resource, pro);
		/********************************************************************************
										   ������Ⱥ�ݻ�
		********************************************************************************/
		//1.�Ӵ���̽���뿪���ı�����2.�Ӵ���̽���뿪���ķ�ʽ
		bool m_evolve_by_predict = false;
		if (!m_evolve_by_predict) {
			bool m_search_balance = true;//ѡ������Ⱥ�Ƿ�E&Eƽ��
			if (!m_search_balance) {
				for (auto& i : pop_exploit_ratio) {
					i = 1.;
				}
			}
			bool m_assign_resource = true;//����Ⱥ�Ƿ������Դ����
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
			//��������Ԥ�⣬����Ԥ�����
			//����Ⱥ�ݻ���ʷ�켣

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
			  ����Ƿ�ϸ���ӿռ�:���оֲ��ṹ���ӿռ䣻λ�ڹȵؽϺõ��ӿռ�
		************************************************************************/
		spaceSubdivision(pro, rnd);

		/**********************************************************************************
							   �ۺϸ�������Ⱥ��Ϣ�������ӿռ���Ϣ
		**********************************************************************************/
		NDSort(offspring_pop);
		SPMOEA::updateVarSpaceInfo(offspring_pop, pro, rnd);
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
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			getPop()[i].envirSelection(pro, rnd, 0);
			//getPop()[i].iteration()++;

			updatePopDistribute(i, pro);
			updatePopdist(i);
		}
		/*********************************************************************************************
			 ���¾��࣬ͳ��һ�����еĵ�ǰ����������ȫΪһ��������С������������Ⱥ��Ȼ������Ⱥ�ݻ�
		**********************************************************************************************/
		updateCluster();

		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		return tag;
	}

	void SPMOEA1_0_1::updatePop(Problem* pro, Algorithm* alg, Random* rnd) {
		auto cur_clusters = getMO_HLC().getClusters();
		//��ȡԭʼ��Ⱥ����
		std::map<size_t, std::vector<std::shared_ptr<Solution<>>>> space_inds;//ÿ���ӿռ京�еĵ�ǰ����
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
		//��ɾ��֮ǰ����Ⱥ
		while (getPop().size() > 0) {
			getPop().remove(getPop().end() - 1);
		}
		//�ٸ��ݾ���Ľ������ԭ�и�����䵽����Ⱥ
		for (size_t i = 0; i < cur_clusters.size(); ++i) {
			size_t count = 0;//�µ����к��еĸ�����
			size_t pop_size = getPopsize();
			for (size_t j = 0; j < cur_clusters[i].size(); ++j) {
				for (size_t k = 0; k < cur_clusters[i][j].size(); ++k) {
					count += space_inds[cur_clusters[i][j][k]].size();
				}
			}
			SPMOEA_pop temp_pop(pop_size, pro);
			if (count > pop_size) {//�����ӿռ������㣬���ӿռ����и�����ѡ��Ӣ����
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
				//�ֲ�������
				auto select_inx = selectIndiFromSpace(temp_pop0, pop_size, pro, rnd);
				for (size_t j = 0; j < select_inx.size(); ++j) {
					temp_pop[j] = temp_pop0[select_inx[j]];
				}
			}
			else {//�����ӿռ䵱ǰ����������ʱ
				size_t cc = 0;//�ӿռ����е�ǰ������
				std::vector<size_t> cluster_space_inx;//�ӿռ����������ӿռ�
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
					//���ӿռ��ڶ������: ��� or ��Ӣ����
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
					//�����ӿռ���Ϣ
					updateVarSpaceInfo(temp_pop1, pro, rnd);
					//������ʷǰ��
					updateHistoryInfo(temp_pop1, pro);
					for (size_t j = 0; j < temp_pop1.size(); ++j) {
						temp_pop[cc + j] = temp_pop1[j];
					}
				}
			}
			//��ʼ���Ӵ�
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
		//��������Ⱥ�ıȽ����������ÿ������Ⱥ��̽���ʺͿ����ʣ��Լ��Ӵ��ĸ���
		std::vector<size_t> pop_min_rank;//����Ⱥ�����rankֵ�Ƚ�
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t min_rank = INT16_MAX;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				if (min_rank > getPop()[i][j].fitness()) {
					min_rank = getPop()[i][j].fitness();
				}
			}
			pop_min_rank.push_back(min_rank);
		}
		std::vector<Real> pop_first_ratio;//����Ⱥǰ�Ÿ���ռ��
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
			Real temp = 1 - 0.2 * pop_min_rank[i];//ϵ��0.5�����ڼ�����Բ��뿪��
			if (temp < 0) {
				temp = 0;
			}
			temp *= 0.5;//ϵ��0.9�������޿����ʣ�0.9*max[0,1-0.5*min_rank]
			pop_exploit_ratio.push_back(temp);
		}
		//�ӿռ串���ʻ�������Ⱥ�ı��־���������Դ���䣬������Ⱥ�����Ӵ��ĸ���
		//����Cָ����㲻ͬ����Ⱥ�����������Ⱥ��֧����
		Matrix dominance_ratio_matrix(getPop().size(), getPop().size());
		for (size_t i = 0; i < dominance_ratio_matrix.row(); ++i) {
			for (size_t j = 0; j < dominance_ratio_matrix.col(); ++j) {
				if (i == j) {
					dominance_ratio_matrix[i][j] = 1;
				}
				else {
					//�����j������Ⱥ�б���i������Ⱥ֧��ĸ�����ռj����Ⱥ���������ı���
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
		//����dominanceȡֵ��Χӳ�������Դ
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
		//	if (exploit_num == 0) {//ʹ����Ⱥ������һ��������̽����
		//		exploit_num++;
		//		explore_num--;
		//	}
		//	if (explore_num == 0) {//ʹ����Ⱥ������һ��������̽����
		//		explore_num++;
		//		exploit_num--;
		//	}
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
		//			//ѡ��ǰ���ӿռ����������ӿռ��ڵĸ��彻��
		//			front_pop.setRate(1, 1. / CAST_CONOP(pro)->numberVariables());
		//			front_pop.setEta(10, 10);
		//			std::vector<size_t> p(2);
		//			do {//���ظ�����
		//				do {
		//					p[0] = front_pop.tournamentSelection(pro, rnd);
		//					p[1] = front_pop.tournamentSelection(pro, rnd);
		//				} while (p[1] == p[0]);
		//				front_pop.crossover(p[0], p[1], front_pop.getOffspring()[j], front_pop.getOffspring()[j + 1], pro, rnd);
		//				front_pop.mutate(front_pop.getOffspring()[j], pro, rnd);
		//				front_pop.mutate(front_pop.getOffspring()[j + 1], pro, rnd);
		//				//�Ӵ�Խ�紦��
		//				repairSol(front_pop.getOffspring()[j].variable().vect(), total_spaces, rnd);
		//				repairSol(front_pop.getOffspring()[j + 1].variable().vect(), total_spaces, rnd);
		//			} while (front_pop.getOffspring()[j].variable().vect() == front_pop[p[0]].variable().vect() || \
		//				front_pop.getOffspring()[j].variable().vect() == front_pop[p[1]].variable().vect() || \
		//				front_pop.getOffspring()[j + 1].variable().vect() == front_pop[p[0]].variable().vect() || \
		//				front_pop.getOffspring()[j + 1].variable().vect() == front_pop[p[1]].variable().vect());

		//			front_pop.getOffspring()[j].setCounter(0);
		//			front_pop.getOffspring()[j + 1].setCounter(0);

		//			//ÿ��ȷ��һ���⣬��2���Ӵ������ѡ��һ��
		//			Real rnd_num = rnd->uniform.next();
		//			if (rnd_num > 0.5) {
		//				getPop()[i].getOffspring()[j].variable() = front_pop.getOffspring()[j + 1].variable();
		//			}
		//		}
		//		else {
		//			std::vector<size_t> p(2);
		//			do {//���ظ�����
		//				do {
		//					p[0] = getPop()[i].tournamentSelection(pro, rnd);
		//					p[1] = getPop()[i].tournamentSelection(pro, rnd);
		//				} while (p[1] == p[0]);
		//				getPop()[i].crossover(p[0], p[1], getPop()[i].getOffspring()[j], getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j], pro, rnd);
		//				getPop()[i].mutate(getPop()[i].getOffspring()[j + 1], pro, rnd);
		//				//�Ӵ�Խ�紦��
		//				repairSol(getPop()[i].getOffspring()[j].variable().vect(), total_spaces, rnd);
		//				repairSol(getPop()[i].getOffspring()[j + 1].variable().vect(), total_spaces, rnd);
		//			} while (getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j].variable().vect() == getPop()[i][p[1]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[0]].variable().vect() || \
		//				getPop()[i].getOffspring()[j + 1].variable().vect() == getPop()[i][p[1]].variable().vect());

		//			getPop()[i].getOffspring()[j].setCounter(0);
		//			getPop()[i].getOffspring()[j + 1].setCounter(0);

		//			//ÿ��ȷ��һ���⣬��2���Ӵ������ѡ��һ��
		//			Real rnd_num = rnd->uniform.next();
		//			if (rnd_num > 0.5) {
		//				getPop()[i].getOffspring()[j].variable() = getPop()[i].getOffspring()[j + 1].variable();
		//			}
		//		}

		//	}
		//	//������Ⱥ�����������ڵ��ӿռ��ڽ���̽��
		//	//����̽�������Ǽ����ӿռ�ĸ�����
		//	bool random_sample = false;
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
		//		//	//����͵ĸ����ʵ��ӿռ�
		//		//	auto sp_inx = std::distance(space_cover_rate.begin(),std::min_element(space_cover_rate.begin(),space_cover_rate.end()));
		//		//	sample_spaces.push_back(total_spaces[sp_inx]);
		//		//	//�����ӿռ�ĸ�����
		//		//	space_cover_rate[sp_inx] = (space_split_num[sp_inx] * space_cover_rate[sp_inx] + 1) / space_split_num[sp_inx];
		//		//	select_count++;
		//		//}
		//	}

		//	//���ӿռ��������������ǻ������и�������½�
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
		//				//�������һ����
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
		//					//�Ӵ�Խ�紦��
		//					repairSol(temp_pop.getOffspring()[0].variable().vect(), bound, rnd);
		//					repairSol(temp_pop.getOffspring()[1].variable().vect(), bound, rnd);
		//				} while (temp_pop.getOffspring()[0].variable().vect() == temp_pop[0].variable().vect() || \
		//					temp_pop.getOffspring()[0].variable().vect() == temp_pop[1].variable().vect() || \
		//					temp_pop.getOffspring()[1].variable().vect() == temp_pop[0].variable().vect() || \
		//					temp_pop.getOffspring()[1].variable().vect() == temp_pop[1].variable().vect());

		//				//ÿ��ȷ��һ���⣬��2���Ӵ������ѡ��һ��
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
		//ϸ�ֱ߽��ӿռ�
		bool split_boundary = 0;
		for (size_t i = 0; i < num_spaces; ++i) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(i);
			auto belong_clusters = getMO_HLC().getSubspaceInfo(i).idx_cluster;
			if (split_boundary) {
				if (belong_clusters.size() > 1) {
					//�����ӿռ��С�����Ƿ�ϸ��
					if (space_volume / total_volume > std::pow(0.1, num_var)) {
						//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
						splitSpace(i, 2 * belong_clusters.size(), 0, 0, true, pro, rnd);
					}
				}
			}
		}
		//ϸ��ȫ��ǰ���ӿռ仹������Ⱥǰ���ӿռ�
		bool split_front = 0;
		auto pre_clusters = getMO_HLC().getClusters();
		if (split_front) {
			for (size_t i = 0; i < num_spaces; ++i) {
				auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(i);
				if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
					if (space_volume / total_volume > std::pow(0.05, num_var)) {
						//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
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
						//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
						int dim = findSplitDim(front_spaces[j], pro);
						auto& space_bound = getMO_HLC().subspaceTree().getBox(front_spaces[j]);
						Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
						splitSpace(front_spaces[j], 2, dim, pos, false, pro, rnd);
					}
				}
			}
		}
		//�ھֲ��ṹ��Ѱ���µľֲ��ṹ
		bool find_new_region = 0;
		//������Ⱥ�ȽϵĽ�������λ��ǰ�ص���Ⱥ���ڵ��ӿռ��Ƿ���пɷ���
		if (find_new_region) {
			bool in_optima_space = true;//��ǰ���ӿռ�Ѱ�Ҿֲ��ṹ
			if (in_optima_space) {
				for (size_t i = 0; i < pre_clusters.size(); ++i) {
					if (getMO_HLC().getSubspaceInfo(pre_clusters[i][0][0]).m_best_rank == 0) {
						//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
						for (size_t j = 0; j < pre_clusters[i][0].size(); ++j) {
							auto split_flag = subspaceSeperable(pre_clusters[i][0][j], pro);//�����ӿռ��Ƿ�ɷ�
							if (split_flag > 0) {
								splitSpace(pre_clusters[i][0][j], 8, 0, 0, true, pro, rnd);//ϸ��2����
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
						auto split_flag = subspaceSeperable(pre_clusters[i][0][j], pro);//�����ӿռ��Ƿ�ɷ�
						if (split_flag > 0) {
							splitSpace(pre_clusters[i][0][j], 8, 0, 0, true, pro, rnd);//ϸ��2����
						}
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

	void SPMOEA1_0_1::updateCluster() {
		clusterSubspace();
		//����Ƚ�ÿ���ӿռ估������ĺû������־ֲ��ṹ
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		//���ҳ�����ֵ����Χ�ӿռ�õ��ӿռ�
		std::vector<size_t> local_good_space;
		bool strict_nei = false;
		for (size_t i = 0; i < space_num; ++i) {
			size_t cur_rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
			auto neigh_inx = getMO_HLC().getSubspaceInfo(i).m_sub_neighbors;
			std::list<int> neighs;
			if (strict_nei) {
				//����һ���غ��ʵ��ӿռ�����
				for (auto k : neigh_inx) {
					auto center_box = getMO_HLC().subspaceTree().getBox(i);
					auto compara_box = getMO_HLC().subspaceTree().getBox(k);
					size_t dim = INT16_MAX;
					for (size_t p = 0; p < center_box.size(); ++p) {
						auto bound1 = center_box[p];
						auto bound2 = compara_box[p];
						//�жϸ�ά�Ƿ��ص�
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
		//�Ծֲ�ǰ���ӿռ����
		std::vector<size_t> flag_idx(space_num, 0);
		std::vector<std::vector<size_t>> cluster_idx;
		cluster_idx = clusterFrontSpace(local_good_space);
		for (size_t i = 0; i < cluster_idx.size(); ++i) {
			for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
				flag_idx[cluster_idx[i][j]] = 1;
			}
		}
		//���浱ǰÿ���������������rankֵ
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
		//����ÿ����������rankֵ�Ĵ�С�����ĸ���������ɢһ��
		while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
			//���ҳ�������С���rankֵ���������
			auto itr = std::min_element(cluster_cur_ranks.begin(), cluster_cur_ranks.end());
			auto cluster_inx = std::distance(cluster_cur_ranks.begin(), itr);
			//��ǰ��������ӿռ�����
			std::vector<size_t> current_cluster_inx;
			for (size_t i = 0; i < final_clusters[cluster_inx].size(); ++i) {
				for (size_t j = 0; j < final_clusters[cluster_inx][i].size(); ++j) {
					current_cluster_inx.push_back(final_clusters[cluster_inx][i][j]);
				}
			}
			//Ȼ���ҳ��������������
			auto cur_cluster = final_clusters[cluster_inx].back();//��ǰ�Ƚϲ�
			std::vector<size_t> neigh_inx;//��ǰ�Ƚϲ���������
			for (size_t i = 0; i < cur_cluster.size(); ++i) {
				auto temp_nei = getMO_HLC().getSubspaceInfo(cur_cluster[i]).m_sub_neighbors;
				std::list<int> neighs;
				if (strict_nei) {
					//����һ���غ��ʵ��ӿռ�����
					for (auto k : temp_nei) {
						auto center_box = getMO_HLC().subspaceTree().getBox(cur_cluster[i]);
						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
						size_t dim = INT16_MAX;
						for (size_t p = 0; p < center_box.size(); ++p) {
							auto bound1 = center_box[p];
							auto bound2 = compara_box[p];
							//�жϸ�ά�Ƿ��ص�
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
			//���ж����������ӿռ��Ƿ񱻵�ǰ���֧�䣬���ǣ�������࣬�����ǣ�������
			bool permit_overlap = true; //�Ƿ������ص������ӿռ�
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
						//����һ���غ��ʵ��ӿռ�����
						for (auto k : nei_inx) {
							auto center_box = getMO_HLC().subspaceTree().getBox(neigh_inx[i]);
							auto compara_box = getMO_HLC().subspaceTree().getBox(k);
							size_t dim = INT16_MAX;
							for (size_t p = 0; p < center_box.size(); ++p) {
								auto bound1 = center_box[p];
								auto bound2 = compara_box[p];
								//�жϸ�ά�Ƿ��ص�
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
					for (size_t j = 0; j < neigh_cur_compara.size(); ++j) {//�����ӿռ�rankֵ�Ƚ�
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
			//�����µ�ǰ���������rank
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
		//�鿴ÿ����ĵ�һ���Ƿ���������������ǣ���ϲ�
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
		//�����ӿռ�������
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

	//void SPMOEA1_0_1::clusterSubspace() {//����汾
	//	//ÿ�ξ���ֻ����һ�㣬����������ص�
	//	getMO_HLC().getClusters().clear();
	//	auto space_num = getMO_HLC().numSubspace();
	//	std::vector<size_t> flag_idx(space_num, 0);
	//	//���ҳ�����ֵ�ǰ���ӿռ�
	//	while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
	//		//�ȵõ������ӿռ���õ�����ֵ������
	//		int rank = INT16_MAX;//���һ�����ǰ��rankֵ
	//		std::vector<size_t> best_idx;//��ǰ����������rankֵ���ӿռ�����
	//		for (size_t i = 0; i < space_num; ++i) {
	//			if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank < rank) {
	//				rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
	//			}
	//		}
	//		for (size_t i = 0; i < space_num; ++i) {
	//			if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank == rank)
	//				best_idx.push_back(i);
	//		}
	//		//����Щ�ǰ���ӿռ���о���
	//		std::vector<std::vector<size_t>> cluster_idx;
	//		cluster_idx = clusterFrontSpace(best_idx);

	//		for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//			for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
	//				flag_idx[cluster_idx[i][j]] = 1;
	//			}
	//		}

	//		bool gradual_cluster = true;//�Ƿ������ɢ
	//		bool strict_cluster = true;//����ÿ���ӿռ䣬ֻ����δ�������������ֵȫ��С�ڸ��ӿռ䣬�ż��뵽����
	//		bool permit_overlap = true;//�Ƿ����������ӿռ��ص�
	//		if (gradual_cluster) {
	//			auto current = cluster_idx;
	//			while (!current.empty()) {
	//				std::vector<std::vector<size_t>> temp_current;//������ĵ�ǰ�Ƚϵ��ӿռ�
	//				for (size_t i = 0; i < current.size(); ++i) {
	//					std::vector<size_t> temp1;//ÿһ����ĵ�һ������
	//					if (current[i].size() > 1 || current[i][0] != INT16_MAX) {
	//						for (size_t j = 0; j < current[i].size(); ++j) {
	//							size_t temp_rank = getMO_HLC().getSubspaceInfo(current[i][j]).m_best_rank;
	//							auto temp_neighbor = getMO_HLC().getSubspaceInfo(current[i][j]).m_sub_neighbors;
	//							std::vector<size_t> sub;
	//							for (auto& k : temp_neighbor) {
	//								if (permit_overlap) {//�Ƿ����������ӿռ��ص�
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
	//							if (permit_overlap) {//�Ƿ����������ӿռ��ص�
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
	//	//�����ӿռ�������
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
	//	//ÿ�ξ���ֻ����һ�㣬����������ص�
	//	getMO_HLC().getClusters().clear();
	//	auto space_num = getMO_HLC().numSubspace();
	//	//���ҳ�����ֵ����Χ�ӿռ�õ��ӿռ�
	//	std::vector<size_t> local_good_space;
	//	bool strict_nei = false;
	//	for (size_t i = 0; i < space_num; ++i) {
	//		size_t cur_rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
	//		auto neigh_inx = getMO_HLC().getSubspaceInfo(i).m_sub_neighbors;
	//		std::list<int> neighs;
	//		if (strict_nei) {
	//			//����һ���غ��ʵ��ӿռ�����
	//			for (auto k : neigh_inx) {
	//				auto center_box = getMO_HLC().subspaceTree().getBox(i);
	//				auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//				size_t dim = INT16_MAX;
	//				for (size_t p = 0; p < center_box.size(); ++p) {
	//					auto bound1 = center_box[p];
	//					auto bound2 = compara_box[p];
	//					//�жϸ�ά�Ƿ��ص�
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
	//	//�Ծֲ�ǰ���ӿռ����
	//	std::vector<size_t> flag_idx(space_num, 0);
	//	std::vector<std::vector<size_t>> cluster_idx;
	//	cluster_idx = clusterFrontSpace(local_good_space);
	//	for (size_t i = 0; i < cluster_idx.size(); ++i) {
	//		for (size_t j = 0; j < cluster_idx[i].size(); ++j) {
	//			flag_idx[cluster_idx[i][j]] = 1;
	//		}
	//	}
	//	//���浱ǰÿ���������������rankֵ
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
	//	//����ÿ����������rankֵ�Ĵ�С�����ĸ���������ɢһ��
	//	while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
	//		//���ҳ�������С���rankֵ���������
	//		auto itr=std::min_element(cluster_cur_ranks.begin(),cluster_cur_ranks.end());
	//		auto cluster_inx = std::distance(cluster_cur_ranks.begin(),itr);
	//		//��ǰ��������ӿռ�����
	//		std::vector<size_t> current_cluster_inx;
	//		for (size_t i = 0; i < final_clusters[cluster_inx].size(); ++i) {
	//			for (size_t j = 0; j < final_clusters[cluster_inx][i].size(); ++j) {
	//				current_cluster_inx.push_back(final_clusters[cluster_inx][i][j]);
	//			}
	//		}
	//		//Ȼ���ҳ��������������
	//		auto cur_cluster = final_clusters[cluster_inx].back();//��ǰ�Ƚϲ�
	//		std::vector<size_t> neigh_inx;//��ǰ�Ƚϲ���������
	//		for (size_t i = 0; i < cur_cluster.size(); ++i) {
	//			auto temp_nei = getMO_HLC().getSubspaceInfo(cur_cluster[i]).m_sub_neighbors;
	//			std::list<int> neighs;
	//			if (strict_nei) {
	//				//����һ���غ��ʵ��ӿռ�����
	//				for (auto k : temp_nei) {
	//					auto center_box = getMO_HLC().subspaceTree().getBox(cur_cluster[i]);
	//					auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//					size_t dim = INT16_MAX;
	//					for (size_t p = 0; p < center_box.size(); ++p) {
	//						auto bound1 = center_box[p];
	//						auto bound2 = compara_box[p];
	//						//�жϸ�ά�Ƿ��ص�
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
	//		//���ж����������ӿռ��Ƿ񱻵�ǰ���֧�䣬���ǣ�������࣬�����ǣ�������
	//		bool permit_overlap = true; //�Ƿ������ص������ӿռ�
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
	//					//����һ���غ��ʵ��ӿռ�����
	//					for (auto k : nei_inx) {
	//						auto center_box = getMO_HLC().subspaceTree().getBox(neigh_inx[i]);
	//						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
	//						size_t dim = INT16_MAX;
	//						for (size_t p = 0; p < center_box.size(); ++p) {
	//							auto bound1 = center_box[p];
	//							auto bound2 = compara_box[p];
	//							//�жϸ�ά�Ƿ��ص�
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
	//				for (size_t j = 0; j < neigh_cur_compara.size(); ++j) {//�����ӿռ�rankֵ�Ƚ�
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
	//		//�����µ�ǰ���������rank
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
	//	//�鿴ÿ����ĵ�һ���Ƿ���������������ǣ���ϲ�
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
	//	//�����ӿռ�������
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
		//ÿ�ξ���ֻ����һ�㣬����������ص�
		getMO_HLC().getClusters().clear();
		auto space_num = getMO_HLC().numSubspace();
		std::vector<size_t> flag_idx(space_num, 0);
		std::vector<std::vector<std::vector<size_t>>> final_clusters;
		//���ҳ�����ֵ�ǰ���ӿռ�
		while (std::find(flag_idx.begin(), flag_idx.end(), 0) != flag_idx.end()) {
			//�ȵõ������ӿռ���õ�����ֵ������
			int rank = INT16_MAX;//���һ�����ǰ��rankֵ
			std::vector<size_t> best_idx;//��ǰ����������rankֵ���ӿռ�����
			for (size_t i = 0; i < space_num; ++i) {
				if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank < rank) {
					rank = getMO_HLC().getSubspaceInfo(i).m_best_rank;
				}
			}
			for (size_t i = 0; i < space_num; ++i) {
				if (flag_idx[i] == 0 && getMO_HLC().getSubspaceInfo(i).m_best_rank == rank)
					best_idx.push_back(i);
			}
			//����Щ�ǰ���ӿռ���о���
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

			bool strict_neighbors = false;//�Ƿ��ص�����
			bool strict_cluster = false;//����ÿ���ӿռ䣬ֻ����δ�������������ֵȫ��С�ڸ��ӿռ䣬�ż��뵽����
			bool permit_overlap = true;//�Ƿ����������ӿռ��ص�
			auto clustered = cluster_idx;
			size_t add_count = cluster_idx.size();
			std::vector<bool> cluster_flag(cur_clusters.size(), true);
			while (add_count > 0) {
				add_count = 0;
				for (size_t i = 0; i < cur_clusters.size(); ++i) {
					std::vector<size_t> temp1;//ÿһ����������
					if (cluster_flag[i]) {
						for (size_t j = 0; j < cur_clusters[i].back().size(); ++j) {
							size_t temp_rank = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_best_rank;
							auto temp_neighbor = getMO_HLC().getSubspaceInfo(cur_clusters[i].back()[j]).m_sub_neighbors;
							std::vector<size_t> neighbors;
							for (auto& ss : temp_neighbor) {
								neighbors.push_back(ss);
							}
							if (strict_neighbors) {
								//����һ���غ��ʵ��ӿռ�����
								std::vector<size_t> overlap_neighbor;
								for (auto& k : temp_neighbor) {
									auto center_box = getMO_HLC().subspaceTree().getBox(cur_clusters[i].back()[j]);
									auto compara_box = getMO_HLC().subspaceTree().getBox(k);
									size_t dim = INT16_MAX;
									for (size_t p = 0; p < center_box.size(); ++p) {
										auto bound1 = center_box[p];
										auto bound2 = compara_box[p];
										//�жϸ�ά�Ƿ��ص�
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
								if (permit_overlap) {//�Ƿ����������ӿռ��ص�
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
			//	std::vector<std::vector<size_t>> temp_current;//������ĵ�ǰ�Ƚϵ��ӿռ�
			//	for (size_t i = 0; i < current.size(); ++i) {
			//		std::vector<size_t> temp1;//ÿһ����ĵ�һ������
			//		if (current[i].size() > 1 || current[i][0] != INT16_MAX) {
			//			for (size_t j = 0; j < current[i].size(); ++j) {
			//				size_t temp_rank = getMO_HLC().getSubspaceInfo(current[i][j]).m_best_rank;
			//				auto temp_neighbor = getMO_HLC().getSubspaceInfo(current[i][j]).m_sub_neighbors;
			//				std::vector<size_t> neighbors;
			//				for (auto &ss : temp_neighbor) {
			//					neighbors.push_back(ss);
			//				}
			//				if (strict_neighbors) {
			//					//����һ���غ��ʵ��ӿռ�����
			//					std::vector<size_t> overlap_neighbor;
			//					for (auto& k : temp_neighbor) {
			//						auto center_box = getMO_HLC().subspaceTree().getBox(current[i][j]);
			//						auto compara_box = getMO_HLC().subspaceTree().getBox(k);
			//						size_t dim = INT16_MAX;
			//						for (size_t p = 0; p < center_box.size(); ++p) {
			//							auto bound1 = center_box[p];
			//							auto bound2 = compara_box[p];
			//							//�жϸ�ά�Ƿ��ص�
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
			//					if (permit_overlap) {//�Ƿ����������ӿռ��ص�
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
		//�����ӿռ�������
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

	//���������ӿռ�ľ��ࣺ��������࣬���ǿ���ǰ�ؽ�ķֲ���
	std::vector<std::vector<size_t>> SPMOEA1_0_1::clusterFrontSpace(const std::vector<size_t>& frontspace) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//���
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space = i;
					head_cluster.push_back(frontspace[begin_space]);//�ȼ���һ���ӿռ�
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
								////�����ӿռ���������ж�
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
				if (temp.empty()) {//û���µ��������
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

	void SPMOEA1_0_1::initiObjSpace(Problem* pro) {
		//for (size_t i = 0; i < m_pop->size(); ++i) {
		//	m_pop->getOffspring()[i] = m_pop->at(i);
		//	m_pop->getParent()[i] = m_pop->at(i);
		//	m_pop->getOffspring()[i + m_pop->size()] = m_pop->at(i);
		//}
		///****************************/
		///*      ����Ŀ��ֵ��Χ      */
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
		///*      ������Ŀ������      */
		///****************************/
		//for (size_t i = 0; i < obj_num; ++i) {
		//	m_subobj_opt_sol.emplace_back(m_pop->at(i));
		//}
		//updateSubObjOpt(m_pop->getParent());
		///****************************/
		///*    Ŀ��ռ仮�ֳ�ʼ��    */
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
